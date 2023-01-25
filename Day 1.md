# Day-1
 ## Introduction to linux

 We have learned here how to do the followings:

 * Command line tool 
 * Basic linux commands

To `make a new folder` use this code:

```
mkdir test
```

This is how er make a github readme file. 
----



I need to change something to commit!


# All I learned Day 1

shell = terminal 
pwd = current directory
ls = alle things, that are current in the directory 
cd = change directory

mkdir -name- = New Folder

always use _ and numbers at the end 
rm = remove
mk= make

rm -r = remove directory
clear terminal 
pfeiltaste hoch/runter = alle commands

touch file = file erstellen 
ls -l= wann files erstellt wurden (Details)
ls -a = shows everything 

rm *.txt = remove all txt-files
ls -l *.txt

touch test.txt --> File erstellt

strg l = delete all

mv file.txt file1.txt ( Name ändern auch mit mv)

back from folder to Desktop --> cd 

rm ./test_001/* .fastq ./test_002/* .fastq --> alle Files with fastq gelöscht

cd filename --> geht in den Ordner

mv = move
mv file.txt ./test/ = move file to /foldername/

cp = copy the file
cp -r test_1/ test_2 (Folder test_1 in Folder test_2 eingefügt)

was ist in dem Folder?   ls ./test_002/

cp ./test_001/*.fastq ./test_002/ --> 

cat >> file.txt --> dann kann man im Dokument schreiben

strg C beendet dann

combine both files --> cat test.txt test2.txt

wget ordnername --> erstellung ordner te

check size of file --> ls -l

unzip al file --> gunzip filename (ersten Buchstaben schreiben und dann mit space auswählen)

mit cat lesen

sudo apt-get install -softwarename- --> Update

User 228
open terminal 
ssh -X sunamXXX@caucluster.rz.uni-kiel.de

passwort: biol217_2022

da 2 Drives : Home and Work

workquota --> zeigt, wieviel Platz ich da habe

cd $WORK --> komme so auf WORK

cd bedeutet nur, dass wir wieder auf HOME zurück kehren 
Auf Home wird alles installiert

alle Folder in WORK

data will be provoided in /home/sunam228

connecting to server
Server:caucluster.rz.uni-kiel.de Port:22
Type:SSH
Folder:/home/sunam228    / 

home is the active path

programming language/what I remember 
R. plot ()
sort.data
ggplot()
data.date


# THis code is about linux introduction 
# this will make a folder

open a document in Visual Sutudio Code with Markdowm All in One 
alles, was dann da geschrieben wird, wird in Github angezeigt
# -> Überschrift 
## -> kleine Schrift

`folder`--> wird farbig hervorgehoben

```
mdkdir test
```
--> Dannw wird das in Code geschrieben

Markdowm cheat sheet 

**Bold text**

*italic*




####Github####

github installieren durch sudo-apt-get install git 

gucken, ob es wirklich installiert wurde : git ---version

git init -->


git status --> 

(base) kurs@Kurs015:~/Desktop/biol217/github_repo$ git config --global user.email "LBertelson29"
(base) kurs@Kurs015:~/Desktop/biol217/github_repo$ git config --global user.email "stu212375@mail.uni-kiel.de
> "

git log --> zeigt an, wer was wie und wo geändert hat im Scriptgco

Resporatory in shell einfügen --> git clone https://....

conda enc list --> how much environments are there at the moment? 

conda list --> how much are there at all?

conda deactivate & conda activate

creating new environment --> conda -n name 

von einem env ins andere --> conda activate name

install new software in new env --> conda install r bei google 

env + software --> conda create -n python_env python

conda env list --> alle env

weg von der base --> conda deactivate


ChatGPT AI kann code erklären 
