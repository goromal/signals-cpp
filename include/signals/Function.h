#pragma once
#include "signals/Signal.h"

using namespace Eigen;

template<typename FuncSpec>
struct Function
{
    using InputType  = typename FuncSpec::InputSignalType::BaseType;
    using OutputType = typename FuncSpec::OutputSignalType::BaseType;

    OutputType operator()(const InputType& input) const
    {
        return FuncSpec::Function(input);
    }

    MatrixXd J(const InputType& input) const
    {
        MatrixXd jacobian;
        return jacobian; // TODO implement; need to take chart map into account
        // TODO maybe move this out into its own vector/scalar-only function
    }

    InputType inputCorrespondingTo(const OutputType& output) const
    {
        InputType input;
        return input; // TODO implement
    }
};
