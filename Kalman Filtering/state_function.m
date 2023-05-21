function new_state = state_function(state, input, A, B)
    new_state = A * state + B * input;
end