#This script produces GC/coverage plots. As input, it uses for each metagenome: 
#1.metaspades-produced nucleotide assemblies (scaffolds.fasta)
#2. CONCOCT-produced coverage files (concoct_depth.txt)
#3. CONCOCT-produced binning files (clustering_gt1000_merged.csv)

library(tidyverse)
library(plotly)
library(RColorBrewer)
library(Biostrings)

#create functions
getGC<-function(fasta){
  #process fasta
  contig_id<-names(fasta) #get contig names
  nucl_freq<-alphabetFrequency(fasta)
  nucl_freq<-as.data.frame(nucl_freq) %>% mutate(gc=100*(G+C+S)/(A+T+W+G+C+S)) #get gc content
  dataset1<-data.frame(contig_id,nucl_freq$gc)
  return(dataset1)
}

######process data

##TS1974
#get coverage
depth1<-read.delim("data/TS1974/concoct_depth.txt")
colnames(depth1)[2]<-"coverage"
depth1<-transform(depth1, contig_id = sub("(.*_cov_[0-9]+\\.[0-9]+).*", "\\1", contig))
depth2<-depth1 %>% group_by(contig_id) %>% summarize(mean_contig=mean(coverage))
#get binning
cluster<-read.csv("data/TS1974/clustering_gt1000_merged.csv")
df<-left_join(depth2,cluster)

#get GC%
fasta<-readDNAStringSet('data/TS1974/scaffolds.fasta')
gc<-getGC(fasta)

df2<-left_join(df,gc)

#filter out short contigs <=1000
df2<-transform(df2, length=sub(".*_length_([0-9]+).*", "\\1", contig_id))
df2$length<-as.numeric(as.character(df2$length))
df2 <- df2 %>% filter(length>999)
colnames(df2)<-c("contig_id","mean_contig","cluster_id","nucl_freq.gc","length")

#make gc_cov plot with contigs colored according to their bins
gc_plot_TS1984<-plot_ly(data=df2,x=~nucl_freq.gc,y=~mean_contig,color=~as.factor(cluster_id),
                 marker=list(size=4)) %>%add_markers()%>%  
  layout( xaxis = list(title = 'GC%'), 
          yaxis = list(title = 'Coverage',type="log"))
htmlwidgets::saveWidget(as.widget(gc_plot_TS1984), "/full/path/gc_cov_plots/prelim_plots/TS1974_gc_plot.html", selfcontained = FALSE) 

#make a dummy variable for graph colors (we used the graph above to identify bins constituiting the lecanoromycete MAG (see methods). Also see methods for how we located ATP9 and mitochondrion contigs)
df2$color<-'1other'
df2$color[df2$cluster_id %in% c(8,39,22)]<-"ascomycete"
df2$color[df2$contig_id=="NODE_1106_length_32113_cov_1337.483000"]<-"mitochondrion"
df2$color[df2$contig_id=="NODE_478_length_80156_cov_71.494426"]<-"atp9"
write.table(df2,"prelim_plots/TS1974_graph.txt",sep='\t',quote = F,row.names = F)



##GTX0161
#get coverage
depth1<-read.delim("data/GTX0161/concoct_depth.txt")
colnames(depth1)[2]<-"coverage"
depth1<-transform(depth1, contig_id = sub("(.*_cov_[0-9]+\\.[0-9]+).*", "\\1", contig))
depth2<-depth1 %>% group_by(contig_id) %>% summarize(mean_contig=mean(coverage))
#get binning
cluster<-read.csv("data/GTX0161/clustering_gt1000_merged.csv")
df<-left_join(depth2,cluster)

#get GC%
fasta<-readDNAStringSet('data/GTX0161/scaffolds.fasta')
gc<-getGC(fasta)

df2<-left_join(df,gc)

#filter out short contigs <=1000
df2<-transform(df2, length=sub(".*_length_([0-9]+).*", "\\1", contig_id))
df2$length<-as.numeric(as.character(df2$length))
df2 <- df2 %>% filter(length>999)
colnames(df2)<-c("contig_id","mean_contig","cluster_id","nucl_freq.gc","length")

#make gc_cov plot with contigs colored according to their bins
gc_plot_GTX0161<-plot_ly(data=df2,x=~nucl_freq.gc,y=~mean_contig,color=~as.factor(cluster_id),
                        marker=list(size=4)) %>%add_markers()%>%  
  layout( xaxis = list(title = 'GC%'), 
          yaxis = list(title = 'Coverage',type="log"))
htmlwidgets::saveWidget(as.widget(gc_plot_GTX0161), "/full/path/gc_cov_plots/prelim_plots/GTX0161_gc_plot.html", selfcontained = FALSE) 


#make a dummy variable for graph colors (we used the graph above to identify bins constituiting the lecanoromycete MAG (see methods). Also see methods for how we located ATP9 and mitochondrion contigs)
df2$color<-'1other'
df2$color[df2$cluster_id %in% c(25,31,18)]<-"ascomycete"
df2$color[df2$contig_id=="NODE_511_length_77571_cov_7767.432259"]<-"mitochondrion"
df2$color[df2$contig_id=="NODE_2671_length_6298_cov_296.856319"]<-"atp9"
write.table(df2,"prelim_plots/GTX0161_graph.txt",sep='\t',quote = F,row.names = F)


##GT0163
#get coverage
depth1<-read.delim("data/GTX0163/concoct_depth.txt")
colnames(depth1)[2]<-"coverage"
depth1<-transform(depth1, contig_id = sub("(.*_cov_[0-9]+\\.[0-9]+).*", "\\1", contig))
depth2<-depth1 %>% group_by(contig_id) %>% summarize(mean_contig=mean(coverage))

#get binning
cluster<-read.csv("data/GTX0163/clustering_gt1000_merged.csv")
df<-left_join(depth2,cluster)

#get GC%
fasta<-readDNAStringSet('data/GTX0163/scaffolds.fasta')
gc<-getGC(fasta)

df2<-left_join(df,gc)

#filter out short contigs <=1000
df2<-transform(df2, length=sub(".*_length_([0-9]+).*", "\\1", contig_id))
df2$length<-as.numeric(as.character(df2$length))
df2 <- df2 %>% filter(length>999)
colnames(df2)<-c("contig_id","mean_contig","cluster_id","nucl_freq.gc","length")

#make gc_cov plot with contigs colored according to their bins
gc_plot_GTX0163<-plot_ly(data=df2,x=~nucl_freq.gc,y=~mean_contig,color=~as.factor(cluster_id),
                         marker=list(size=4)) %>%add_markers()%>%  
  layout( xaxis = list(title = 'GC%'), 
          yaxis = list(title = 'Coverage',type="log"))
htmlwidgets::saveWidget(as.widget(gc_plot_GTX0163), "/full/path/gc_cov_plots/prelim_plots/GTX0163_gc_plot.html", selfcontained = FALSE) 


#make a dummy variable for graph colors (we used the graph above to identify bins constituiting the lecanoromycete MAG (see methods). Also see methods for how we located ATP9 and mitochondrion contigs)
df2$color<-'1other'
df2$color[df2$cluster_id %in% c(55,67,19,82)]<-"ascomycete"
df2$color[df2$contig_id=="NODE_339_length_76204_cov_3544.912553"]<-"mitochondrion"
df2$color[df2$contig_id=="NODE_161_length_119058_cov_145.352495"]<-"atp9"
write.table(df2,"prelim_plots/GTX0163_graph.txt",sep='\t',quote = F,row.names = F)

##GTX0158
#get coverage
depth1<-read.delim("data/GTX0158/concoct_depth.txt")
colnames(depth1)[2]<-"coverage"
depth1<-transform(depth1, contig_id = sub("(.*_cov_[0-9]+\\.[0-9]+).*", "\\1", contig))
depth2<-depth1 %>% group_by(contig_id) %>% summarize(mean_contig=mean(coverage))

#get binning
cluster<-read.csv("data/GTX0158/clustering_gt1000_merged.csv")
df<-left_join(depth2,cluster)

#get GC%
fasta<-readDNAStringSet('data/GTX0158/scaffolds.fasta')
gc<-getGC(fasta)

df2<-left_join(df,gc)

#filter out short contigs <=1000
df2<-transform(df2, length=sub(".*_length_([0-9]+).*", "\\1", contig_id))
df2$length<-as.numeric(as.character(df2$length))
df2 <- df2 %>% filter(length>999)
colnames(df2)<-c("contig_id","mean_contig","cluster_id","nucl_freq.gc","length")

#make gc_cov plot with contigs colored according to their bins
gc_plot_GTX0158<-plot_ly(data=df2,x=~nucl_freq.gc,y=~mean_contig,color=~as.factor(cluster_id),
                         marker=list(size=4)) %>%add_markers()%>%  
  layout( xaxis = list(title = 'GC%'), 
          yaxis = list(title = 'Coverage',type="log"))
htmlwidgets::saveWidget(as.widget(gc_plot_GTX0158), "~/Documents/gulya/mito_mystery/git_version/gc_cov_plots/prelim_plots/GTX0158_gc_plot.html", selfcontained = FALSE) 


#make a dummy variable for graph colors (see methods for how we located the bin constituiting the ascomycete MAG and for how we located ATP9 and mitochondrion contigs)
df2$color<-'1other'
df2$color[df2$cluster_id ==70]<-"ascomycete"
df2$color[df2$contig_id=="NODE_254_length_28338_cov_10948.085139"]<-"mitochondrion"
df2$color[df2$contig_id=="NODE_4_length_1304943_cov_350.589746"]<-"atp9"
write.table(df2,"prelim_plots/GTX0158_graph.txt",sep='\t',quote = F,row.names = F)


###Make the plot
#collect the data
TS1974<-read.delim("prelim_plots/TS1974_graph.txt")
GTX0158<-read.delim("prelim_plots/GTX0158_graph.txt")
GTX0161<-read.delim("prelim_plots/GTX0161_graph.txt")
GTX0163<-read.delim("prelim_plots/GTX0163_graph.txt")
TS1974$metagenome<-"TS1974"
GTX0158$metagenome<-"GTX0158"
GTX0161$metagenome<-"GTX0161"
GTX0163$metagenome<-"GTX0163"
plot_data<-rbind(GTX0158,GTX0161,TS1974,GTX0163)
plot_data$metagenome<-factor(plot_data$metagenome,levels=c("TS1974","GTX0163","GTX0158","GTX0161"))


#make FINAL faceted plot
palette<-c('#b5b8b5','#ea9e0f','#7a03bd','#ee041b')

gc_plot<-plot_data%>%
  group_by(metagenome) %>%
  do(p=plot_ly(., x = ~nucl_freq.gc, y = ~mean_contig, color = ~as.factor(color),colors = palette, type = "scatter",marker=list(size=3)) %>%
       layout( xaxis = list(title = 'GC%',range = c(0,100),tickvals = list(0, 25, 50,75,100)), 
               yaxis = list(title = 'Coverage',type="log"))) %>% 
  subplot(nrows = 2, shareX = TRUE, shareY = TRUE)

htmlwidgets::saveWidget(as.widget(gc_plot), "/full/path/gc_cov_plots/final_plot/gc_plot.html", selfcontained = FALSE) 

Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/opt/anaconda2/bin/", sep = .Platform$path.sep))
orca(gc_plot,"final_plot/gc_plot.svg")

