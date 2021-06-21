# kallistobustools

kallisto | bustools workflow for pre-processing single-cell RNA-seq data

Please visit https://kallistobus.tools for tutorials on how to process single-cell RNA-seq data.

## Contributing

Create a Google Colab notebook and make a pull request. If approved `your_tutorial.ipynb` will be added to a new folder `tutorials/docs/tutorials/your_tutorial_folder/{python,R}`. 

You will also need to make a pull request to add a reference to your tutorial in the `mkdocs.yml` file. For example:

```yaml
- Tutorials: 
  - Your Tutorial Name:
    - Python: tutorials/your_tutorial_folder/{python,R}/your_tutorial.md
```

*Note*: tutorials will either be in python or R, so paths should be written as `tutorials/your_tutorial_folder/R/your_tutorial.md` for example.
