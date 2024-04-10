library("hexSticker")

`quad.form` <- function (M, x){drop(crossprod(crossprod(M,Conj(x)),x))}



rot <- function(theta){
  matrix(c(
    cos(theta),-sin(theta),
    sin(theta), cos(theta)
  ),2,2,byrow=TRUE)
}

squash <- function(alpha,beta){diag(c(alpha,beta))}

ellipses <- function(alpha,beta,theta){
    quad.form(squash(alpha,beta),rot(theta))
}


plotter <- function(...){
  par(mar=c(5,2,4,2)+0) 
  par(pty="s")
  M <- ellipses(2,0.5,pi/6)
  x <- seq(from=-1,to=1,by=0.01)
  y <- x
  p <- as.matrix(expand.grid(x,y))
  v <- apply(p,1,function(x){quad.form(M,x)})
  contour(x,y,matrix(v,length(x),length(y)),
        levels=c(0.01,0.1,0.3,0.6,1),
        lwd=16,drawlabels=FALSE)
}

plotter()
bmp(file="quadform_icon.bmp",bg="#7733FF",width=2000,height=1500)
plotter()
dev.off()

sticker("quadform_icon.bmp", package="quadform", p_size=28, s_x=1, s_y=0.99,
        s_width=2,asp=sqrt(3)/2, white_around_sticker=TRUE, h_fill="#7733FF",
        h_color="#000000", filename="quadform.png")

