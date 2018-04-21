
library(dplyr, warn = FALSE)

# Must manually set
multiplier = 5

stat_file = file.path("..", "hetnets", "permuted", "stats.tsv")
stat_df = readr::read_tsv(stat_file) %>%
  dplyr::mutate(complete = round(complete * multiplier, 2))
head(stat_df, 4)

unchanged_df = stat_df %>%
  # Average over permutations
  dplyr::group_by(abbrev, complete) %>%
  dplyr::summarize(unchanged = mean(unchanged)) %>%
  dplyr::ungroup() %>%
  dplyr::bind_rows(dplyr::data_frame(abbrev=unique(stat_df$abbrev), complete = 0, unchanged = 1))

abbrevs = unchanged_df %>%
  dplyr::filter(complete == multiplier) %>%
  dplyr::arrange(desc(unchanged)) %>%
  .[['abbrev']]

unchanged_df %>%
  ggplot2::ggplot(ggplot2::aes(x = complete, y = 100 * unchanged, color = abbrev)) +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::scale_colour_discrete(breaks = abbrevs, name='Metaedge') +
  ggplot2::xlab('Attempt multiplier') +
  ggplot2::ylab('Percent of Edges Unchanged')

bar_df = stat_df %>%
  tidyr::gather(key = 'measure', value = 'percent', duplicate:excluded, same_edge:undirected_duplicate) %>%
  dplyr::group_by(abbrev, measure) %>%
  dplyr::summarize(
    percent = 100 * weighted.mean(percent, attempts)
  ) %>%
  dplyr::filter(measure != 'excluded')

bar_df$abbrev = factor(bar_df$abbrev, levels=abbrevs)

bar_df %>%
  dplyr::filter(measure %in% c('duplicate')) %>%
  ggplot2::ggplot(ggplot2::aes(x = abbrev, y = percent, fill = measure)) +
  ggplot2::geom_bar(stat = "identity", position = "dodge") +
  ggplot2::coord_flip()

bar_df %>%
  dplyr::filter(!(measure %in% c('duplicate', 'unchanged'))) %>%
  ggplot2::ggplot(ggplot2::aes(x = abbrev, y = percent, fill=measure)) +
  ggplot2::geom_bar(stat = "identity", position = "dodge") +
  ggplot2::coord_flip()
