---
layout: page
title: "Find the secret BUS file"
---

{% include JB/setup %}

The image you clicked on has a __secret__! Follow the instructions below to find out more:

### 1. Download the image
Right-click the tsne on the previous page and click "Copy Link Address". Navigate to your terminal and download the image using
```
$ wget https://www.kallistobus.tools/assets/secret_tsne.jpg
```
### 2. Install the unrar utility
If you have homebrew installed on your machine then run 
```
$ brew install unrar
```
If not, then [install homebrew](https://docs.brew.sh/Installation) and run the command above.

### 3. Un-zip the image (you read that right!)
Run 
```
$ unzip secret_tsne.jpg
```

### 5. Un-rar the rar file
With the unrar utility, run
```
$ unrar e secret.rar 
```

to find the surprise :).