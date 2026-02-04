# FEM Theory

Theory of Finite Element Method

## MathJax Test Page

When \(a \ne 0\), there are two solutions to \(ax^2 + bx + c = 0\) and they are

\begin{equation}
ax^2 + bx + c = 0
\end{equation}

$$ x = {-b \pm \sqrt{b^2-4ac} \over 2a} $$

## state diagrams with mermaid

``` mermaid
stateDiagram-v2
  state fork_state <<fork>>
    [*] --> fork_state
    fork_state --> State2
    fork_state --> State3

    state join_state <<join>>
    State2 --> join_state
    State3 --> join_state
    join_state --> State4
    State4 --> [*]
```

## Using class diagrams

``` mermaid
classDiagram
  Person <|-- Student
  Person <|-- Professor
  Person : +String name
  Person : +String phoneNumber
  Person : +String emailAddress
  Person: +purchaseParkingPass()
  Address "1" <-- "0..1" Person:lives at
  class Student{
    +int studentNumber
    +int averageMark
    +isEligibleToEnrol()
    +getSeminarsTaken()
  }
  class Professor{
    +int salary
  }
  class Address{
    +String street
    +String city
    +String state
    +int postalCode
    +String country
    -validate()
    +outputAsLabel()
  }
```