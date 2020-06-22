# Matrix

A Matrix helper library(?)
I was coding some neural networks in c++ and i need a faster way to make matrixes and operations so i decided to create some helper functions and where we are.

I plan to upgrade this library and make actually a util thing

## Current state

### Declaration

It can't be more simple:

    Matrix m; // Empty matrix 0x0 no data
    Matrix m(3x3); // Matrix 3x3 with 0 as initial value
    Matrix m(3,3,5);  // Matrix 3x3 with 5 as initial value

### Operations

- Scalar

        m += 10; // Adding 10 in each element
        m = m * 20; // Mutiply by 20 each element

- Other matrix

        m += m; // Adding two Matrixes
        m -= m; // Complex statemens are not recomended
        m *= m; // In special when dealing with mutiplications

## Misc

### AI Realated

- Randomize
  Generating random values using a Mersenne Twister

        m.randomize(); // Values of m are random values, in range [-1,1]

- Map

        m.map(fun); // Apply an function in every element of m

        double fun(double x){
                reutrn exp(x);
        }

## Issues

Found something? Any bug? See space to improvement? Make a issue and/or a pull request, fell free to help.

I'm very appreciated for any help

Made with ❤ and C++