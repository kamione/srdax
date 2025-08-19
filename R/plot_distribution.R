plot_null_distribution <- function(dat) {

    if (dat$p_value == 0) {
        dat$p_value = "< 0.001"
    } else {
        dat$p_value = paste0("= ", round(dat$p_value, 3))
    }

    annotation_text <- paste0(
        "*β* = ",
        round(dat$empirical_rho, 3),
        ", *p* ",
        dat$p_value
    )

    fig <- data.frame(x = dat$permuted_rhos) %>%
        ggplot(aes(x = x)) +
        geom_density(color = "gray70", fill = "gray70") +
        geom_vline(
            xintercept = dat$empirical_rho,
            color = "tomato3",
            linetype = "dashed",
            linewidth = 1
        ) +
        labs(
            x = "Standardized Coefficients (*β*)",
            y = "",
            title = annotation_text
        ) +
        ggthemes::theme_pander() +
        theme(
            plot.margin = margin(5, 5, 5, 5, "mm"),
            plot.title = ggtext::element_markdown(),
            axis.title.x = ggtext::element_markdown(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
    return(fig)
}
