---
title: Home
layout: home
nav_order: 1
---
![](https://github.com/Dragonmasterx87/CompBio1-Tulane/blob/main/assets/images/logo_small.png)
## Welcome to the Computational Biology Workshop 1 (CBW1) workspace
### Change website Color scheme <button class="btn js-toggle-dark-mode">Change Website to Dark Color</button>

<script>
const toggleDarkMode = document.querySelector('.js-toggle-dark-mode');

jtd.addEvent(toggleDarkMode, 'click', function(){
  if (jtd.getTheme() === 'dark') {
    jtd.setTheme('light');
    toggleDarkMode.textContent = 'Change to Darth Vader';
  } else {
    jtd.setTheme('dark');
    toggleDarkMode.textContent = 'Change to Anakin Skywalker';
  }
});
</script>

### Introduction
In this workshop attendees will learn the basics of scRNAseq experimental design. We will be using the Seurat package in R to analyze single cell RNA sequencing data.
This workshop will outline the basics of R syntax and usage using R and the integrated developmental environment RStudio. This course is designed to be a basic introduction,
where the fundamentals of scRNAseq analysis will be outlined. 

### Prerequisites
_Nulla scientia necessaria_. It is expected that you have no prior coding eperience and want to learn how to analyze scRNAseq data from scratch, if you have coding experience then
this course may act as a refresher to go over the basic workflow of Seurat. It is advised to bring your own laptop with a functional operating system so that you may follow the instructor's
workflow.

### Learning Outline
### Day 1
1. Introduction to single cell RNA sequencing
2. Experimental design
3. Comparison of library preparation and sequencing platforms
4. Overview of single-cell analysis using Seurat
5. Installing R and basic syntax
 
### Day 2
1. Data loading/QC/filtering
2. Normalization using VST and SC Transform
3. Multidimensional reduction and cell clustering
4. Differential testing framework
5. Visualization
6. Pitfalls, caveats, and future topics

### Project
In order to be eligible for the attendance certificate you must complete the project, dont worry at the end of this workshop you will be able to perform the necessry computation to complete this.


----
[Just the Docs]: https://just-the-docs.github.io/just-the-docs/
[GitHub Pages]: https://docs.github.com/en/pages
[README]: https://github.com/just-the-docs/just-the-docs-template/blob/main/README.md
[Jekyll]: https://jekyllrb.com
[GitHub Pages / Actions workflow]: https://github.blog/changelog/2022-07-27-github-pages-custom-github-actions-workflows-beta/
[use this template]: https://github.com/just-the-docs/just-the-docs-template/generate
