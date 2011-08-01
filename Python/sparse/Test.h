#include <iostream>


namespace test {

class SomeClass {
public:
    SomeClass() {
        std::cout << "TestClass constructor called" << std::endl;
    }

    virtual ~SomeClass() {
        std::cout << "TestClass destructor called" << std::endl;
    }

    int method1() {
        std::cout << "TestClass::method1 called" << std::endl;
        return 23;
    }
};

}
