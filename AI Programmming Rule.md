### **Programming Guidelines**

Please adhere to the following programming guidelines:

1.  **Project Goal**: To build a general-purpose 3D Finite Element Method (FEM) simulator for coupled electro-thermal problems.

2.  **Core Design Principles (Open/Closed Principle)**:
    *   **Extension Method**: Add new functionality by creating new derived classes and implementing the **existing** virtual functions of the base class.
    *   **Interface Stability**: The public interface of base classes is considered stable and **must not** be modified to support new features. If a new function must be added to the base class for logical reasons, it must be explicitly identified.
    *   **Instance Creation**: Create instances of derived classes via a Factory Pattern, which returns a smart pointer to the base class (e.g., `std::unique_ptr<Base>`).

3.  **Interface Rules**:
    *   Avoid adding new `public` functions to derived classes.
    *   If a new function is absolutely necessary, it must be explicitly marked with an `// [EXTENSION]` comment.

4.  **Commenting Standard**:
    *   Use Doxygen-style comments.
    *   When modifying code, you must update the corresponding comments to keep them in sync.
    *   
5.  **Code Review** All code changes must undergo a code review process to ensure adherence to these guidelines.
   
6.  **Interaction Principle**: use English when you are coding，and Chinese for communication and discussion with user, as well as commentings.