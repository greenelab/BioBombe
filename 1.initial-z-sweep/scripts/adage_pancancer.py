"""
Gregory Way 2017
Variational Autoencoder - Pan Cancer
scripts/adage_pancancer.py

Comparing a VAE learned features to ADAGE features. Use this script within
the context of a parameter sweep to compare performance across a grid of
hyper parameters.

Usage:

    Run in command line with required command arguments:

        python scripts/adage_pancancer.py --learning_rate
                                          --batch_size
                                          --epochs
                                          --sparsity
                                          --noise
                                          --output_filename
                                          --num_components
                                          --scale
                                          --subset_mad_genes
                                          --dataset

    Typically, arguments to this script are compiled automatically.

    See `scripts/num_components_paramsweep.py` for more details

Output:
Loss and validation loss for the specific model trained
"""

import os
import argparse
import numpy as np
import pandas as pd

from sklearn.preprocessing import MinMaxScaler
from keras.engine.topology import Layer
from keras.layers import Input, Dense, Dropout, Activation
from keras.models import Sequential, Model
import keras.backend as K
from keras.regularizers import l1
from keras import optimizers, activations


class TiedWeightsDecoder(Layer):
    """
    Transpose the encoder weights to apply decoding of compressed latent space
    """
    def __init__(self, output_dim, encoder, activation=None, **kwargs):
        self.output_dim = output_dim
        self.encoder = encoder
        self.activation = activations.get(activation)
        super(TiedWeightsDecoder, self).__init__(**kwargs)

    def build(self, input_shape):
        self.kernel = self.encoder.weights
        super(TiedWeightsDecoder, self).build(input_shape)

    def call(self, x):
        # Encoder weights: [weight_matrix, bias_term]
        output = K.dot(x - self.encoder.weights[1],
                       K.transpose(self.encoder.weights[0]))
        if self.activation is not None:
            output = self.activation(output)
        return output

    def compute_output_shape(self, input_shape):
        return (input_shape[0], self.output_dim)


parser = argparse.ArgumentParser()
parser.add_argument('-l', '--learning_rate',
                    help='learning rate of the optimizer')
parser.add_argument('-b', '--batch_size',
                    help='Number of samples to include in each learning batch')
parser.add_argument('-e', '--epochs',
                    help='How many times to cycle through the full dataset')
parser.add_argument('-s', '--sparsity',
                    help='How much L1 regularization penalty to apply')
parser.add_argument('-n', '--noise',
                    help='How much Gaussian noise to add during training')
parser.add_argument('-f', '--output_filename',
                    help='The name of the file to store results')
parser.add_argument('-c', '--num_components', default=100,
                    help='The latent space dimensionality to test')
parser.add_argument('-o', '--optimizer', default='adam',
                    help='optimizer to use', choices=['adam', 'adadelta'])
parser.add_argument('-w', '--untied_weights', action='store_false',
                    help='use tied weights in training ADAGE model')
parser.add_argument('-a', '--scale', action='store_true',
                    help='Add decision to scale input data')
parser.add_argument('-m', '--subset_mad_genes', default=8000,
                    help='The number of mad genes to subset')
parser.add_argument('-x', '--dataset', default='TCGA',
                    choices=['TCGA', 'TARGET', 'GTEX'],
                    help='the dataset to use in the sweep')
args = parser.parse_args()

# Set hyper parameters
learning_rate = float(args.learning_rate)
batch_size = int(args.batch_size)
epochs = int(args.epochs)
sparsity = float(args.sparsity)
noise = float(args.noise)
output_filename = args.output_filename
latent_dim = int(args.num_components)
use_optimizer = args.optimizer
tied_weights = args.untied_weights
scale = args.scale
subset_mad_genes = args.subset_mad_genes
dataset = args.dataset

# Random seed
seed = int(np.random.randint(low=0, high=10000, size=1))
np.random.seed(seed)

# Load Data
file = 'train_{}_expression_matrix_processed.tsv.gz'.format(dataset.lower())
rnaseq_file = os.path.join('..', '0.expression-download', 'data', file)
rnaseq_df = pd.read_table(rnaseq_file, index_col=0)

# Determine most variably expressed genes and subset
if subset_mad_genes is not None:
    mad_genes = rnaseq_df.mad(axis=0).sort_values(ascending=False)
    top_mad_genes = mad_genes.iloc[0:subset_mad_genes, ].index
    rnaseq_df = rnaseq_df.loc[:, top_mad_genes]

# Zero One normalize input data
if scale:
    scaler = MinMaxScaler()
    x = scaler.fit_transform(rnaseq_df)
    rnaseq_df = pd.DataFrame(x, index=rnaseq_df.index,
                             columns=rnaseq_df.columns)

original_dim = rnaseq_df.shape[1]

# Split 10% test set randomly
test_set_percent = 0.1
rnaseq_test_df = rnaseq_df.sample(frac=test_set_percent)
rnaseq_train_df = rnaseq_df.drop(rnaseq_test_df.index)

if tied_weights:
    # Input place holder for RNAseq data with specific input size
    encoded_rnaseq = Dense(latent_dim,
                           input_shape=(original_dim, ),
                           activity_regularizer=l1(sparsity),
                           activation='relu')
    dropout_layer = Dropout(noise)
    tied_decoder = TiedWeightsDecoder(input_shape=(latent_dim, ),
                                      output_dim=original_dim,
                                      activation='sigmoid',
                                      encoder=encoded_rnaseq)

    autoencoder = Sequential()
    autoencoder.add(encoded_rnaseq)
    autoencoder.add(dropout_layer)
    autoencoder.add(tied_decoder)

else:
    input_rnaseq = Input(shape=(original_dim, ))
    encoded_rnaseq = Dropout(noise)(input_rnaseq)
    encoded_rnaseq_2 = Dense(latent_dim,
                             activity_regularizer=l1(sparsity))(encoded_rnaseq)
    activation = Activation('relu')(encoded_rnaseq_2)
    decoded_rnaseq = Dense(original_dim, activation='sigmoid')(activation)
    autoencoder = Model(input_rnaseq, decoded_rnaseq)

if use_optimizer == 'adadelta':
    optim = optimizers.Adadelta(lr=learning_rate)
elif use_optimizer == 'adam':
    optim = optimizers.Adam(lr=learning_rate)

autoencoder.compile(optimizer=optim, loss='mse')

hist = autoencoder.fit(np.array(rnaseq_train_df), np.array(rnaseq_train_df),
                       shuffle=True,
                       epochs=epochs,
                       batch_size=batch_size,
                       validation_data=(np.array(rnaseq_test_df),
                                        np.array(rnaseq_test_df)))

# Save training performance
history_df = pd.DataFrame(hist.history)
history_df = history_df.assign(num_components=latent_dim,
                               learning_rate=learning_rate,
                               batch_size=batch_size,
                               epochs=epochs,
                               sparsity=sparsity,
                               noise=noise,
                               seed=seed,
                               dataset=dataset)
history_df.to_csv(output_filename, sep='\t')
