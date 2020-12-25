## to be sourced in results.R


ggp_theme_default <- theme(
    panel.background=element_rect(fill="white"),
    panel.grid.major=element_line(color="white"),
    panel.grid.minor=element_line(color="white"),
    plot.background=element_rect(fill="white"),
    plot.margin=margin(t=1, r=1, b=1, l=1, unit="lines"),
    plot.title=element_text(size=12, face="bold", hjust=0.5),
    plot.subtitle=element_text(size=10, face="plain", hjust=0.5),
    axis.title.x=element_text(size=10, face="plain"),
    axis.title.y=element_text(size=10, face="plain"),
    axis.text.x=element_text(size=9, face="plain"),
    axis.text.y=element_text(size=9, face="plain"),
    # axis.ticks.x=element_none(),
    # axis.ticks.y=element_none(),
    axis.line.x.bottom=element_line(),
    axis.line.y.left=element_line(),
    legend.key=element_rect(fill="white"),
    legend.position="top",
    legend.title=element_blank(),
    strip.background=element_rect(fill="black"),
    strip.text=element_text(colour="white")
)


ggp_theme_col <- ggp_theme_default +
    theme(
        plot.margin=margin(t=0.5, r=0.5, b=0.5, l=0.5, unit="lines"),
        plot.title=element_text(size=12, face="plain", hjust=0.5),
        axis.text.x=element_text(size=10, face="plain", angle=45, hjust=1),
        axis.text.y=element_text(size=10, face="plain"),
        axis.line.x.bottom=element_blank(),
        axis.line.y.left=element_blank(),
        legend.position="none"
    )


ggp_theme_vol <- ggp_theme_default +
    theme(
        plot.margin=margin(t=0.25, r=0.25, b=0.25, l=0.25, unit="lines"),
        legend.key=element_blank(),
        legend.key.size=unit(0.001, "cm"),
        legend.key.width=unit(0.001, "cm"),
        legend.position=c(0.18, 0.9),
        legend.spacing.x=unit(0.001, "cm"),
        legend.spacing.y=unit(0.001, "cm"),
        legend.text=element_text(size=8, margin=margin(t=0.001))
    )


ggp_theme_box <- ggp_theme_default +
    theme(
        plot.margin=margin(t=0.5, r=1, b=0.5, l=1, unit="lines"),
        axis.text.x=element_text(size=9, face="plain", angle=45, hjust=1),
        axis.text.y=element_text(size=10, face="plain"),
        legend.position="none",
        legend.key=element_blank(),
        legend.key.size=unit(0.75, "cm"),
        legend.key.width=unit(0.75, "cm"),
        legend.spacing.x=unit(0.5, "cm"),
        legend.spacing.y=unit(0.5, "cm"),
        legend.text=element_text(size=9, margin=margin(t=0.1))
    )


draw_key_default <- function(data, params, size) {
    ## this function is here in case GeomBar$draw_key must be manipulated and
    ## then reset once occasion for altered parameters done
    ## https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2
    if (is.null(data$size)) {
        data$size <- 0.5
    }
    lwd <- min(data$size, min(size) / 4)
    grid::rectGrob(
        width = grid::unit(1, "npc") - grid::unit(lwd, "mm"),
        height = grid::unit(1, "npc") - grid::unit(lwd, "mm"),
        gp = grid::gpar(
            col = data$colour %||% NA,
            fill = alpha(data$fill %||% "grey20", data$alpha),
            lty = data$linetype %||% 1,
            lwd = lwd * .pt,
            linejoin = params$linejoin %||% "mitre",
            lineend = if (identical(params$linejoin, "round")) "round" else "square"
        )
    )
}


draw_key_large <- function(data, params, size) {
    ## this function is here in case GeomBar$draw_key must be reset to allow for
    ## larger spacing between legend keys; not used in this script
    ## https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2
    lwd <- min(data$size, min(size) / 4)
    grid::rectGrob(
        width = grid::unit(0.6, "npc"),
        height = grid::unit(0.6, "npc"),
        gp = grid::gpar(
            col = data$colour,
            fill = alpha(data$fill, data$alpha),
            lty = data$linetype,
            lwd = lwd * .pt,
            linejoin = "mitre"
        )
    )
}
