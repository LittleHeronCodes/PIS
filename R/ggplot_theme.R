#` ggplot Transparent theme
#`
#` Theme for transparent ggplot figures. function adopts from theme_grey.
#' @param base_size font base size. default 11
#' @export

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

