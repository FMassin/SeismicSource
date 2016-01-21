# SeismicSource
A simple module for basic seismic source modeling.

## Importing the SeismicSource class
Start python and run:
  ```python
  from source import SeismicSource
  ```
The `SeismicSources` class then be used.

## Example
Import the SeismicSource class and run:
  ```python
  model_example = SeismicSource([1,2,3,4,5,6])
  model_example.Aki_Richards.plot('P')
  ``` 

## General installation 
Start python and run:
  ```python
  import source
  ```
The whole module is then accessible by typing `source.` and TAB.
