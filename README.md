<!-- Your Project title, make it sound catchy! -->

# Navier-Stokes Propagator

<!-- Provide a short description to your project -->

## Description

The repository contains source files for a Julia package written to solve the chaotic Navier-Stokes equations. The project aims to exemplify the improvement of simulation code quality due to the use of good programming abstractions (Julia Structs), done by defining sensible simulation constructs that reflects each core component of the simulation. Additionally, the repository adopts good programming practice through implementing detailed documentation with tests written for each simulation construct.

The algorithm here used is a replicate of Diablo (https://github.com/johnryantaylor/DIABLO), written by John R. Taylor, albeit with a minor modification to the mesh spacing used. 

<!-- What should the students going through your exemplar learn -->

## Learning Outcomes

- Ability to distinguish and compartmentalise a multi-facet simulation algorithm into sensible simulation constructs via the implementation with Julia Structs 
- Adopting good testing methodologies, through implementation of unit and integration tests
- Understanding of the FFTW and HDF5 library and their abstraction through the usage of utility modules in their respective source files

<!-- How long should they spend reading and practising using your Code.
Provide your best estimate -->

| Task       | Time    |
| ---------- | ------- |
| Reading    | 3 hours |
| Practising | 3 hours |

## Requirements

<!--
If your exemplar requires students to have a background knowledge of something
especially this is the place to mention that.

List any resources you would recommend to get the students started.

If there is an existing exemplar in the ReCoDE repositories link to that.
-->

### Theoretical
- Understanding of the Fourier series and the relation to the discrete Fourier transform
- Understanding of the finite difference method for approximating spatial derivatives and the use of implicit / explicit time-stepping schemes
- Familiarity with Julia, readers not familiar with Julia can refer the ECRI Julia course (https://www.imperial.ac.uk/students/academic-support/graduate-school/professional-development/doctoral-students/research-computing-data-science/courses/introduction-to-julia/)

Note : The theoratical understanding of the finite difference method, Fourier series and the time-stepping scheme are not a major necessity. An understanding of the mathematical structure for algorithimic implementation is sufficient to understand the ideas discussed in this tutorial.
<!-- List the system requirements and how to obtain them, that can be as simple
as adding a hyperlink to as detailed as writting step-by-step instructions.
How detailed the instructions should be will vary on a case-by-case basis.

Here are some examples:

- 50 GB of disk space to hold Dataset X
- Anaconda
- Python 3.11 or newer
- Access to the HPC
- PETSc v3.16
- gfortran compiler
- Paraview
-->

### System

- Julia 1.11 with HDF5, FFTW and Test packages installed
<!-- Instructions on how the student should start going through the exemplar.

Structure this section as you see fit but try to be clear, concise and accurate
when writing your instructions.

For example:
Start by watching the introduction video,
then study Jupyter notebooks 1-3 in the `intro` folder
and attempt to complete exercise 1a and 1b.

Once done, start going through through the PDF in the `main` folder.
By the end of it you should be able to solve exercises 2 to 4.

A final exercise can be found in the `final` folder.

Solutions to the above can be found in `solutions`.
-->

## Getting Started

1. Clone the repository
2. Launch Julia pkg mode via "]" in the REPL.
3. Activate the package environment via "activate ." .
4. Testing can be done by entering "pkg> test" .
<!-- An overview of the files and folder in the exemplar.
Not all files and directories need to be listed, just the important
sections of your project, like the learning material, the code, the tests, etc.

A good starting point is using the command `tree` in a terminal(Unix),
copying its output and then removing the unimportant parts.

You can use ellipsis (...) to suggest that there are more files or folders
in a tree node.

-->

## Project Structure

```log
.
├── examples
│   ├── ex1
│   └── ex2
├── src
|   ├── file1.py
|   ├── file2.cpp
|   ├── ...
│   └── data
├── app
├── docs
├── main
└── test
```

<!-- Change this to your License. Make sure you have added the file on GitHub -->

## License

This project is licensed under the [BSD-3-Clause license](LICENSE.md)
