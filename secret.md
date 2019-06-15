---
layout: page
title: "Find the BUS"
---

{% include JB/setup %}

The image you clicked on has a __secret__! Follow the instructions below to find it:

### 1. Download the image
Right-click the image on the previous page and download the image. Navigate to the terminal and download the image using
```
$ wget https://www.kallistobus.tools/assets/secret_tsne.jpg
```
### 2. Install the unrar utility
If you have homebrew installed on your machine then run 
```
$ brew install unrar
```
If not, then [install homebrew](https://brew.sh/) and run the command above. 

### 3. Unzip the image
Run 
```
$ unzip secret_tsne.jpg
```

### 5. Unrar the rar file
With the unrar utility, run
```
$ unrar e secret.rar 
```

### 6. Gunzip the file
```
$ gunzip output.sort.bus
```
### 7. Find the BUS!
```
$ bustools text -p output.sort.bus | head
```
