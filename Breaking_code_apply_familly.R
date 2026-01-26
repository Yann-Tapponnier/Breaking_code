############################################################################
#### Apply familly ####
############################################################################


Sure! The **`apply` family** refers to a group of functions in **R** (and conceptually in other languages like Python too) used for **applying a function** to elements of a data structure (like vectors, lists, matrices, or data frames) **without using explicit loops**.


### Summary Table

| Function   | Input            | Output        | Key Use Case                            |
  | ---------- | ---------------- | ------------- | --------------------------------------- |
  | `apply()`  | matrix/array     | vector/matrix | Row/column operations                   |
  | `lapply()` | list/vector      | list          | Apply function, keep as list            |
  | `sapply()` | list/vector      | simplified    | Apply function, simplify result         |
  | `tapply()` | vector + factor  | array         | Grouped summary stats                   |
  | `mapply()` | multiple vectors | vector/list   | Apply function element-wise in parallel |
  
  
  Here’s a brief overview of the main `apply` family in **R**:
  
  ---
  
  ### 📦 `apply()`
  
  * **Used for**: matrices or arrays
* **Purpose**: Apply a function to rows or columns
* **Example**:
  
  
  apply(matrix, 1, sum)  # sum of each row
apply(matrix, 2, mean) # mean of each column


---
  
  ### 🧩 `lapply()`
  
  * **Used for**: lists or vectors
* **Returns**: a **list**
  * **Example**:
  
  
  lapply(list(1:5, 6:10), mean)


---
  
  ### 🧮 `sapply()`
  
  * **Used for**: lists or vectors
* **Returns**: a **simplified result** (vector, matrix) if possible
* **Example**:
  
  
  sapply(list(1:5, 6:10), mean)


---
  
  ### 📘 `tapply()`
  
  * **Used for**: applying a function over **groups**
  * **Example**:
  
  
  tapply(iris$Sepal.Length, iris$Species, mean)


---
  
  ### 🗂️ `mapply()`
  
  * **Used for**: applying a function to **multiple arguments in parallel**
  * **Example**:
  
  
  mapply(sum, 1:3, 4:6)  # sum(1,4), sum(2,5), sum(3,6)