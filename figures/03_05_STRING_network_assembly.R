# STRING network assembly

library(rsvg)
library(grImport2)
library(grid)
library(gridExtra)
library(xml2)
library(stringr)

px_to_pt <- function(px) px * 72/96
pt_to_px <- function(pt) pt * 96/72

bump_svg_typography <- function(in_svg, out_svg, min_pt = 8, stroke_scale = 1) {
  doc  <- read_xml(in_svg)
  root <- xml_root(doc)

  # ensure a viewBox exists
  if (is.na(xml_attr(root, "viewBox"))) {
    w <- xml_attr(root, "width"); h <- xml_attr(root, "height")
    if (!is.na(w) && !is.na(h) && grepl("px$", w) && grepl("px$", h)) {
      vw <- as.numeric(sub("px$", "", w)); vh <- as.numeric(sub("px$", "", h))
      xml_set_attr(root, "viewBox", sprintf("0 0 %g %g", vw, vh))
    }
  }

  min_px <- pt_to_px(min_pt)

  # --- helper to bump a single numeric value with unit
  bump_val <- function(val_chr, min_pt, min_px) {
    m <- str_match(val_chr, "^\\s*([0-9.]+)\\s*(pt|px)?\\s*$")
    if (any(is.na(m))) return(val_chr)
    val  <- as.numeric(m[2]); unit <- m[3]
    if (!is.na(unit) && unit == "pt") {
      if (val < min_pt) sprintf("%gpt", min_pt) else sprintf("%gpt", val)
    } else {
      # px or unitless -> treat as px
      if (val < min_px) sprintf("%gpx", min_px) else sprintf("%gpx", val)
    }
  }

  # --- 1) bump font-size in style="" blocks
  style_nodes <- xml_find_all(doc, "//*[@style]")
  for (n in style_nodes) {
    style <- xml_attr(n, "style")
    kv <- unlist(strsplit(style, ";", fixed = TRUE))
    kv <- kv[nzchar(trimws(kv))]
    if (length(kv)) {
      parts <- strsplit(kv, ":", fixed = TRUE)
      keys  <- trimws(vapply(parts, `[`, "", 1))
      vals  <- trimws(vapply(parts, function(x) if (length(x) >= 2) paste(x[-1], collapse=":") else "", ""))
      # bump font-size if present
      if (any(keys == "font-size")) {
        idx <- which(keys == "font-size")
        vals[idx] <- vapply(vals[idx], bump_val, character(1), min_pt = min_pt, min_px = min_px)
      }
      # (optional) scale stroke-width
      if (!isTRUE(all.equal(stroke_scale, 1)) && any(keys == "stroke-width")) {
        idx <- which(keys == "stroke-width")
        vals[idx] <- vapply(vals[idx], function(v) {
          m <- str_match(v, "^\\s*([0-9.]+)\\s*(pt|px)?\\s*$")
          if (any(is.na(m))) return(v)
          newv <- as.numeric(m[2]) * stroke_scale
          unit <- m[3]; if (!nzchar(unit)) unit <- "px"
          sprintf("%g%s", newv, unit)
        }, character(1))
      }
      # reassemble
      new_style <- paste0(keys, ": ", vals, collapse = "; ")
      xml_set_attr(n, "style", new_style)
    }
  }

  # --- 2) bump font-size attributes on <text> / <tspan>
  textish <- xml_find_all(doc, "//svg:text|//svg:tspan", xml_ns(doc))
  for (n in textish) {
    fs <- xml_attr(n, "font-size")
    if (!is.na(fs)) {
      xml_set_attr(n, "font-size", bump_val(fs, min_pt, min_px))
    }
  }

  # --- 3) also scale stroke-width attributes if requested
  if (!isTRUE(all.equal(stroke_scale, 1))) {
    stroke_nodes <- xml_find_all(doc, "//*[@stroke-width]")
    for (n in stroke_nodes) {
      sw <- xml_attr(n, "stroke-width")
      m <- str_match(sw, "^\\s*([0-9.]+)\\s*(pt|px)?\\s*$")
      if (!any(is.na(m))) {
        newv <- as.numeric(m[2]) * stroke_scale
        unit <- m[3]; if (!nzchar(unit)) unit <- "px"
        xml_set_attr(n, "stroke-width", sprintf("%g%s", newv, unit))
      }
    }
  }

  write_xml(doc, out_svg)
  invisible(out_svg)
}

in1 <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig4_STRING/CD4.svg"
in2 <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig4_STRING/CD8.svg"
clean1 <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig4_STRING/CD4_clean.svg"
clean2 <- "/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig4_STRING/CD8_clean.svg"

# 1) Normalize the SVGs (librsvg simplifies CSS etc.)
rsvg::rsvg_svg(in1, clean1)
rsvg::rsvg_svg(in2, clean2)

# 2) Import as 'Picture' (now parseable)
big1 <- sub("\\.svg$", "_min6pt.svg", clean1)
big2 <- sub("\\.svg$", "_min6pt.svg", clean2)
bump_svg_typography(clean1, big1, min_pt = 6, stroke_scale = 1)
bump_svg_typography(clean2, big2, min_pt = 6, stroke_scale = 1)

# import the *bumped* SVGs (not the originals)
pic1 <- grImport2::readPicture(big1)
pic2 <- grImport2::readPicture(big2)

g1 <- grImport2::pictureGrob(pic1)
g2 <- grImport2::pictureGrob(pic2)

# 3) Lay out two networks on a 6" x 7" page (vector output)
# target canvas
W <- 7; H <- 6
# choose exact widths so both fit (leave a 0.2" gutter)
left_w  <- unit(3.45, "in")
gap_w   <- unit(0.20, "in")
right_w <- unit(3.45, "in")

svglite::svglite("/Chlamydia-SingleCell/SCTransformV2 Oct26 2023/Figures/Fig4_STRING/two_networks.svg", width = W, height = H)
grid.newpage()

# 3-column layout: left figure, gutter, right figure
pushViewport(viewport(
  layout = grid.layout(
    nrow = 1, ncol = 3,
    widths  = unit.c(left_w, gap_w, right_w),
    heights = unit(H, "in")
  )
))

## LEFT cell — clip ON; draw pic1 at explicit width in inches
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, clip = "on"))
grid.picture(pic1, x = 0.5, y = 0.5, just = "center",
             width = left_w)               # <- exact physical width
popViewport()

## RIGHT cell
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3, clip = "on"))
grid.picture(pic2, x = 0.5, y = 0.5, just = "center",
             width = right_w)
popViewport(2)
dev.off()

