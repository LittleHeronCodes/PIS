## ggplot theme for transparent figures

#' Transparent theme for ggplot
#'
#' Theme for transparent ggplot figures. function adopts from theme_grey.
#' @param base_size base font size. default 11
#' @param base_family base font family.
#' @param base_line_size base line size. default base_size / 22
#' @param base_rect_size base rectangle size. default base_size / 22
#' @return ggplot theme
#' @import ggplot2
#' @export
#' @examples
#' ggplot2::theme_set(theme_transparent3)

theme_transparent3 <- function(base_size = 11, base_family = "", base_line_size = base_size/22, base_rect_size = base_size/22) {
	
	.theme <- theme_grey(base_size = base_size, base_family = base_family, 
		base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
	theme(
		# theme_bw
		panel.border = element_rect(fill = NA, colour = "grey20"), 
		panel.grid = element_line(colour = "grey92"), 
		panel.grid.minor = element_line(size = rel(0.5)), 

		# transparent background
		plot.background = element_rect(fill = "transparent",colour = NA),
		panel.background = element_rect(fill = "transparent", colour = NA),
		strip.background = element_rect(fill = "transparent", colour = NA, size = 0.7), 
		legend.background = element_rect(fill = "transparent", colour = NA), 
		legend.key = element_blank(), 
		panel.ontop = FALSE, complete = TRUE
		)
}

