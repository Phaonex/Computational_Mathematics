module probability
using DataFrames

@doc """
    event_occurency(event, Es) -> Vector | Number

Compute the frequecy of the `event` ∊ occurency given an `Es` ∊ Ω space.


# Examples
```julia-repl
julia> event_occurency(3,  rand(1:6, 10))
1
```
"""
event_occurency(event::Union{Number, String}, Es:: Union{Tuple, Vector, Array}) =
    count(e -> e == event , Es)

"""
    prob(Es::Union{Tuple, Vector, Array}, Ω=1) -> Vector | Number

Compute the probability of the `Es` ∊ Ω in `Ω` ∊ Ω space.


# Examples
```julia-repl
julia> event_occurency(3,  rand(1:6, 10))
1
```
"""    
prob(Es::Union{Number, Tuple, Vector, Array}, Ω=1) = Es / Ω

const face_cards = ['J', 'K', 'Q', 'A']
const suits = ['♣', '♢', '♡', '♠']
const numbers = [2:10;]
const colors = ["black", "red"]
const ranks = union(numbers, face_cards)

const datas = (
    Coin = Set(['H', 'T']), # assume that heads = 1 and tails = 0. = Ω space
    Dice = Set(1:6;), # Ω space
    cards_df = DataFrame(suits = suits, color=repeat(colors, 2), ranks=map(c -> map(s -> s, ranks), suits)),
    
)
end # module probability
