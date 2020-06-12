//Automatic templates instanciation with variadic templates metaprogramming
//cf https://stackoverflow.com/questions/25202250/c-template-instantiation-avoiding-long-switches#
//Default case:
template <int ...> struct IntList {};
int handle_cases(int dimensions, int argc, const char ** argv, IntList<>)
{
    std::cerr << "The " << dimensions << " dimensions case is not staticly implemented." << std::endl;
    std::cerr << "Add it to the template in Tools/dimensionsInstanciation.hpp at line 33." << std::endl;
    exit(1);
}
//Recursive pass on variadic template parameter:
template <int I, int ...N> int handle_cases(int dimensions, int argc, const char ** argv, IntList<I, N...>)
{
    if(I != dimensions) return handle_cases(dimensions, argc, argv, IntList<N...>());
    return main_template<VecX<I>>(argc, argv);
}
template <int ...N> int handle_cases(int dimensions, int argc, const char ** argv)
{
    return handle_cases(dimensions, argc, argv, IntList<N...>());
}


int main(int argc, const char ** argv)
{
    int dimensions = DIM;
    for(int i(0); i < argc; i++)
        if(strcmp(argv[i], "-d") == 0 && i + 1 < argc)
        {
            dimensions = std::stoi(argv[i + 1]);
            break;
        }
    return handle_cases<2, 3, 4, 5, 6, 7, 8, 20>(dimensions, argc, argv);
}
