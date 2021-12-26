com_data=read.delim("S:/2018/pglemasle/Data/data_com.txt")


levels(com_data$Milieux)


Mare=com_data[com_data[,1]=="Mare",]

Mare=droplevels(Mare)


levels(Mare$Id_Sites)



Mare[,2]=as.character(Mare[,2], trim = TRUE)

#Mare[,2]=as.numeric(Mare[,2])

levels(as.factor(Mare[,2]))


#Mare[,2]=paste("M", Mare[,2], sep="_")

Mare[,2]

#gsub(" ", "0", Mare[,2], fixed = TRUE)

Mare[,2]

Mare[,2]=as.factor(Mare[,2])


#####

envi_data=read.delim("S:/2018/pglemasle/Data/envi_data.txt", header=TRUE)

envi_data[,2]=gsub("M_", "", envi_data[,2], fixed = TRUE)

envi_data[,2]

envi_data[,2]=sub("^0+", "", envi_data[,2])

envi_data[,2]

####

test=as.character(Mare[!duplicated(Mare[,2]),2])

!(test %in% envi_data[,2])

test[!(test %in% envi_data[,2])]

length(test[!(test %in% envi_data[,2])])

####

Mare_select=Mare[(test %in% envi_data[,2]),]


####

write.table(Mare_select, "Mare_select.txt" , sep="\t",dec=".")

write.table(Mare_select, "Mare_select.txt" , sep="\t",dec=".")

####


