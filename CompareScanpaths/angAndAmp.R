
theta = hist(180*saccDat$dir/pi,45)
saccAngle = data.frame(theta=theta$mids%%360, count=theta$density)
saccAmp = data.frame(amp=saccDat$amp/60)

ggplot(saccAmp, aes(x=amp))+ geom_histogram() + scale_x_continuous(name="saccadic amplitude", limits=c(0,20))+ scale_y_continuous(name="number of saccades")+theme_bw()
ggsave("saccAmp.pdf", width=5, height=5)



roseplt <- ggplot(saccAngle, aes(x=theta, y=count)) + geom_bar(width=10, stat="identity") 
roseplt <- roseplt + scale_x_continuous(name=" ", limits = c(0, 360),breaks = c(0, 45, 90, 135, 180, 225, 270, 315))+scale_y_continuous(name=" ",breaks=NULL)+ coord_polar(direction=-1)+theme_bw()
ggsave("roseplot.pdf", width=5, height=5)


