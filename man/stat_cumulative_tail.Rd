% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_cumulative_tails.R
\name{stat_cumulative_tail}
\alias{stat_cumulative_tail}
\title{Plot Cumulative Tail}
\usage{
stat_cumulative_tail(
  mapping = NULL,
  data = NULL,
  geom = "point",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  bins = 100,
  ...
)
}
\arguments{
\item{mapping}{Set of aesthetic mappings created by \code{\link[ggplot2:aes]{aes()}}. If specified and
\code{inherit.aes = TRUE} (the default), it is combined with the default mapping
at the top level of the plot. You must supply \code{mapping} if there is no plot
mapping.}

\item{data}{The data to be displayed in this layer. There are three
options:

If \code{NULL}, the default, the data is inherited from the plot
data as specified in the call to \code{\link[ggplot2:ggplot]{ggplot()}}.

A \code{data.frame}, or other object, will override the plot
data. All objects will be fortified to produce a data frame. See
\code{\link[ggplot2:fortify]{fortify()}} for which variables will be created.

A \code{function} will be called with a single argument,
the plot data. The return value must be a \code{data.frame}, and
will be used as the layer data. A \code{function} can be created
from a \code{formula} (e.g. \code{~ head(.x, 10)}).}

\item{geom}{geom}

\item{position}{Position adjustment, either as a string naming the adjustment
(e.g. \code{"jitter"} to use \code{position_jitter}), or the result of a call to a
position adjustment function. Use the latter if you need to change the
settings of the adjustment.}

\item{na.rm}{If \code{FALSE}, the default, missing values are removed with
a warning. If \code{TRUE}, missing values are silently removed.}

\item{show.legend}{logical. Should this layer be included in the legends?
\code{NA}, the default, includes if any aesthetics are mapped.
\code{FALSE} never includes, and \code{TRUE} always includes.
It can also be a named logical vector to finely select the aesthetics to
display.}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetics,
rather than combining with them. This is most useful for helper functions
that define both data and aesthetics and shouldn't inherit behaviour from
the default plot specification, e.g. \code{\link[ggplot2:borders]{borders()}}.}

\item{bins}{number of bins}

\item{...}{Other arguments passed on to \code{\link[ggplot2:layer]{layer()}}. These are
often aesthetics, used to set an aesthetic to a fixed value, like
\code{colour = "red"} or \code{size = 3}. They may also be parameters
to the paired geom/stat.}
}
\description{
Nice extention of the default ggplot2 stats. However, it should be easier to
use the plot_cumulative_tails() function.
}
\details{
\code{stat_bin()} is suitable only for continuous x data. If your x data is
discrete, you probably want to use \code{\link[ggplot2:stat_count]{stat_count()}}.

By default, the underlying computation (\code{stat_bin()}) uses 30 bins;
this is not a good default, but the idea is to get you experimenting with
different number of bins. You can also experiment modifying the \code{binwidth} with
\code{center} or \code{boundary} arguments. \code{binwidth} overrides \code{bins} so you should do
one change at a time. You may need to look at a few options to uncover
the full story behind your data.

In addition to \code{geom_histogram()}, you can create a histogram plot by using
\code{scale_x_binned()} with \code{\link[ggplot2:geom_bar]{geom_bar()}}. This method by default plots tick marks
in between each bar.
}
\section{Orientation}{


This geom treats each axis differently and, thus, can thus have two orientations. Often the orientation is easy to deduce from a combination of the given mappings and the types of positional scales in use. Thus, ggplot2 will by default try to guess which orientation the layer should have. Under rare circumstances, the orientation is ambiguous and guessing may fail. In that case the orientation can be specified directly using the \code{orientation} parameter, which can be either \code{"x"} or \code{"y"}. The value gives the axis that the geom should run along, \code{"x"} being the default orientation you would expect for the geom.

}

\section{Aesthetics}{


\code{geom_histogram()} uses the same aesthetics as \code{\link[ggplot2:geom_bar]{geom_bar()}};
\code{geom_freqpoly()} uses the same aesthetics as \code{\link[ggplot2:geom_line]{geom_line()}}.

}

\section{Computed variables}{


These are calculated by the 'stat' part of layers and can be accessed with \link[ggplot2:aes_eval]{delayed evaluation}.
\itemize{
\item \code{after_stat(count)}\cr number of points in bin.
\item \code{after_stat(density)}\cr density of points in bin, scaled to integrate to 1.
\item \code{after_stat(ncount)}\cr count, scaled to a maximum of 1.
\item \code{after_stat(ndensity)}\cr density, scaled to a maximum of 1.
\item \code{after_stat(width)}\cr widths of bins.
}

}

\section{Dropped variables}{


\describe{
\item{\code{weight}}{After binning, weights of individual data points (if supplied) are no longer available.}
}

}

\examples{
library(ggplot2)
data("tcga_brca_test")
snvs <- SNVs(tcga_brca_test)

ggplot(snvs, aes(VAF, color = sample_id)) +
  stat_cumulative_tail()

ggplot(snvs, aes(VAF, y = stat(y), color = sample_id)) +
  stat_cumulative_tail()
}
\seealso{
\code{\link[ggplot2:stat_count]{stat_count()}}, which counts the number of cases at each x
position, without binning. It is suitable for both discrete and continuous
x data, whereas \code{stat_bin()} is suitable only for continuous x data.
}
