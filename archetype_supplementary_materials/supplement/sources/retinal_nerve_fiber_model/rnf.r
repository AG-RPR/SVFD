# Copyright (C) 2012-2014 Tobias Elze
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


##########
# aux. functions for implementation of the RNF model described by
# Jansonius et al. (2009, 2012):
##########
b.nasal <- function(phi0) 0.00083*phi0**2 + 0.02*phi0 - 2.65
b.nasal.ci.lower <- function(phi0) 0.0005*phi0**2 + 0.011*phi0 - 6.8
b.nasal.ci.upper <- function(phi0) 0.0012*phi0**2 + 0.05*phi0 - 1.51

asup <- function(phi0) (phi0-121)/14
c.superior <- function(phi0) 1.9 + 1.4*tanh(asup(phi0))
b.superior <- function(phi0, betas = -1.9) exp(betas + 3.9*tanh(-asup(phi0)))
b.superior.ci.upper <- function(phi0) b.superior(phi0, betas = -1.3)
b.superior.ci.lower <- function(phi0) b.superior(phi0, betas = -2.5)

ainf <- function(phi0) (-phi0-90)/25
c.inferior <- function(phi0) 1 + 0.5*tanh(ainf(phi0))
b.inferior <- function(phi0, betai = 0.7) -exp(betai+1.5*tanh(-ainf(phi0)))
b.inferior.ci.upper <- function(phi0) b.inferior(phi0, betai = 1.3)
b.inferior.ci.lower <- function(phi0) b.inferior(phi0, betai = 0.1)


#' Calculate rnf object according to Jasonius et al. (2012), Eq. 1
#' from radius and phi0
#'
#' @param r radius
#' @param phi0 angle at the starting point at a circle with radius r0 = 4
#' @param confidence.interval whether to add confidence interval of phi
#'   if TRUE, ci.lower, ci.upper are added
#' @references
#'  Jansonius, M, Schiefer, J., Nevalainen, J., Paetzold, J., Schiefer, U.
#'  A mathematical model for describing the retinal nerve fiber bundle trajectories in
#'  the human eye: Average course, variability, and influence of refraction, optic disc
#'  size and optic disc position
#'  Experimental Eye Research 105 (2012) 70-78
#' @export
rnf <- function(r, phi0, confidence.interval=FALSE, r0=4)
{
	phi0 = ifelse(phi0<=180, phi0, phi0-360)
	fiber = list()
	fiber$phi0 = phi0
	fiber$r = r
	calc.phi <- function(bsup=b.superior, binf = b.inferior, bnas = b.nasal)
	{
		if(phi0>=0)
		{
			par.c = c.superior(phi0)
			par.b = ifelse(
				60 <= phi0, 
				bsup(phi0), 
				bnas(phi0))
		}
		else
		{
			par.c = c.inferior(phi0)
			par.b = ifelse(
				phi0 <= -60, 
				binf(phi0), 
				bnas(phi0))
		}
		phi0 + par.b * (r-r0)**par.c
	}
	fiber$phi = calc.phi()
	if(confidence.interval)
	{
		fiber$ci.lower = calc.phi(b.superior.ci.lower, b.inferior.ci.lower, b.nasal.ci.lower)
		fiber$ci.upper = calc.phi(b.superior.ci.upper, b.inferior.ci.upper, b.nasal.ci.upper)
	}
	attr(fiber, 'class') <- 'rnf'
	fiber
}

#' Tries to estimate the angle phi0 at the starting point at a circle with radius r0 = 4 
#' from a given radius r and angle phi; inverse function of Jasonius et al. (2012), Eq. 1. 
#' Note: r and phi are given in Jansonius' coordinate sytem, not in coordinate system
#' of the visual field (see \code{\link{phi0.at.vf.location}} for the latter).
#'
#' @param r radius
#' @param phi angle
#' @seealso \code{\link{rnf}}, \code{\link{phi0.at.vf.location}}
#' @export
rnf.phi0 <- function(phi, r)
{
	phi = ifelse(phi<=180, phi, phi-360)
	iv = if (phi>=0) c(0, 180) else c(-180, 0)
	uniroot(function(x) rnf(r, phi0=x)$phi - phi, iv)$root
}

#' Tries to estimate the angle phi0 (see \code{\link{rnf}} for details) for a location
#' in the visual field, given either in polar (r, phi) or cartesian coordinates (x, y),
#' from a given radius r and angle phi; inverse function of Jasonius et al. (2012), Eq. 1. 
#' Note: the nerve fiber coordinate system is mirrored on the x axis relative to the VF
#'
#' @param r radius (for polar coordinates, ignored if x and y are given)
#' @param phi angle (for polar coordinates, ignored if x and y are given)
#' @param x VF x location in cartesian coordinates; causes to ignore phi and r
#' @param y VF y location in cartesian coordinates; causes to ignore phi and r
#' @seealso \code{\link{rnf}}, \code{\link{rnf.phi0}}
#' @export
phi0.at.vf.location <- function(phi=45, r=sqrt(18), x=NULL, y=NULL)
{
	if(is.null(x))
	{
		x = r*cos(pi*phi/180)
		y = r*sin(pi*phi/180)
	}
	x.jans = x - 15
	# the nerve fiber coordinate system is mirrored on the x axis relative to the VF:
	y.jans = -y - 2
	r.jans = sqrt(x.jans**2 + y.jans**2)
	phi.jans = 180*atan2(y.jans, x.jans)/pi
	rnf.phi0(phi.jans, r.jans)
}


#' Convert polar coordinate system centered on optic disc
#' to cartesian system centered on fovea
#'
#' @param r radius
#' @param phi angle
#' @references
#'  Jansonius, M, Schiefer, J., Nevalainen, J., Paetzold, J., Schiefer, U.
#'  A mathematical model for describing the retinal nerve fiber bundle trajectories in
#'  the human eye: Average course, variability, and influence of refraction, optic disc
#'  size and optic disc position
#'  Experimental Eye Research 105 (2012) 70-78
#' @export
rnf2cartesian <- function(r, phi) 
	list(
		x = r * cos(pi*phi/180) + 15,
		y = r * sin(pi*phi/180) + 2) 

#' plot retinal nerve fibres given in polar coordinates centered on optic disc
#' on a coordinate system centered on the fovea;
#' requires library plotrix
#'
#' @param rnf rnf object (as returned by function rnf)
#' @param confidence.interval whether to add the confidence interval
#'   (only applies if phi is a list as returned by \code{\link{rnf.phi}})
#' @param sector which sector ("nasal", "inferior", or "superior") the RNF
#'   is supposed to belong to; by default guessed from rnf$phi0
#' @param col color of the RNF line
#' @param lty line type of the RNF line
#' @param add whether to add to an existing plot
#' @param project.on.vf whether to project RNF trajectory to visual field, i.e.
#'   mirror the plot along the x-axis
#' @param cutat30degree whether to cut the RNF trajectory at 30 degree
#' @param clip.to.sector whether to restrict the RNF trajectory to the sector
#'   to which it is supposed to belong to (option sector)
#' @param annotate.angle angle to annotate on the RNF trajectory;
#'   set NULL to suppress annotation
#' @param xaxis.annot.cex cex of angle labels on x-axis
#' @param ... additional arguments are passed to \code{\link{plot}}
#' @references
#'  Jansonius, M, Schiefer, J., Nevalainen, J., Paetzold, J., Schiefer, U.
#'  A mathematical model for describing the retinal nerve fiber bundle trajectories in
#'  the human eye: Average course, variability, and influence of refraction, optic disc
#'  size and optic disc position
#'  Experimental Eye Research 105 (2012) 70-78
#' @seealso \code{\link{rnf.phi}}
#' @export
plot.rnf <- function(
	rnf, 
	confidence.interval = FALSE,
	sector = NULL, 
	col="blue", 
	lty="solid",
	add=FALSE, 
	project.on.vf = FALSE,
	cutat30degree = TRUE, 
	clip.to.sector = TRUE,
	annotate.angle = rnf$phi0,
	xaxis.annot.cex = 1,
	...)
{
	require("plotrix")
	phi = rnf$phi
	r = rnf$r
	if(is.null(sector)) sector = ifelse(abs(rnf$phi0)<60, "nasal", ifelse(rnf$phi0<0, "inferior", "superior"))
	cart <- rnf2cartesian(r, phi)
	if(!add)
	{
		plot(
			c(0,0), c(-30, 30), 
			col="gray", 
			type="l", 
			xlim=c(-30, 30), ylim=c(-30, 30), 
			xaxt="n", yaxt = "n", 
			xlab="", ylab="",
			bty="n", 
			...)
		lines(c(-30, 30), c(0, 0), col="gray")
		lapply(5*(1:6), function(r) draw.circle(0, 0, r, border="gray"))
		draw.circle(15, ifelse(project.on.vf, -1, 1)*2, 4, col="lightgray", lty="dotted")
	}
	if(cutat30degree)
	{
		within30deg = sqrt(cart$x**2 + cart$y**2) <= 30 & phi <= 360
		cart$x = cart$x[within30deg]
		cart$y = cart$y[within30deg]
		r = r[within30deg]
	}
	if(clip.to.sector && sector != "nasal")
	{
		tolerance.radius = 7
		if(sector=="superior")
			withinsector = cart$y >= 0
		else
			withinsector = r<tolerance.radius | cart$y < 0
		cart$x = cart$x[withinsector]
		cart$y = cart$y[withinsector]
	}
	lines(cart$x, ifelse(project.on.vf, -1, 1)*cart$y, col=col, lty=lty)
	if(confidence.interval)
	{
		upper = rnf
		lower = rnf
		upper$phi=rnf$ci.upper
		lower$phi=rnf$ci.lower
		lapply(
			list(lower, upper), 
			function(v) plot.rnf(
				v,
				sector=sector, 
				col=col, 
				lty="dashed", 
				add=TRUE, 
				project.on.vf = project.on.vf,
				cutat30degree=cutat30degree, 
				clip.to.sector=clip.to.sector,
				annotate.angle = NULL,
				...))
	}
	if(!is.null(annotate.angle))
	{
		annot.x = cart$x[round(length(cart$x)/2)]
		annot.y = cart$y[round(length(cart$x)/2)]
		boxed.labels(
			annot.x, ifelse(project.on.vf, -1, 1)*annot.y, 
			labels=as.expression(bquote(.(annotate.angle)^o)), 
			cex=0.5,
			border=F, 
			col=col)
	}
	lapply(10*(1:3), function(i) boxed.labels(-i, 0, border=F, labels=as.expression(bquote(.(i)^o)), cex=xaxis.annot.cex))
	if(!add) text(15, 2, "ONH", cex=xaxis.annot.cex)
	boxed.labels(0, 0, border=F, "F", cex=xaxis.annot.cex)
	NULL
}


