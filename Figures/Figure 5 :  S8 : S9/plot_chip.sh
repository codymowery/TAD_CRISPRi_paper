## Plotting c(CD28_start-10e3,TAD_stop)
## hg38 chr2:203696475-203967032
## hg19 chr2:204561198-204831755
pyGenomeTracks --tracks chipAtac.ini --region chr2:204561198-204831755 --width 10 --height 5 --trackLabelHAlign center -o ./ChIP_ATAC.pdf

##hg19 liftover of hg38 chr2:203814752-203816578
pyGenomeTracks --tracks zoom_hg19.ini --region chr2:204679475-204681301 --width 10 --height 5 --trackLabelHAlign center -o ./ChIP_Zoom_hg19.pdf
pyGenomeTracks --tracks zoom_hg38.ini --region chr2:203814752-203816578 --width 10 --height 5 --trackLabelHAlign center -o ./ChIP_Zoom_hg38.pdf
