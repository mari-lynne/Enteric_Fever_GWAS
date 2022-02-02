#SET PLINK TO PATH ####
#Save Plink in folder in PATH directory so can be called dircetly

echo $PATH
#add line to shell start up file
export PATH=$PATH:/place/with/the/file

#my version
#edit the .bashrc file in vim and add PATH file line in Vim
vim ~/.bashrc
export PATH="$PATH:/home/mari/Plink"
source ~/.bashrc