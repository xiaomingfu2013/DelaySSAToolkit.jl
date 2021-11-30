"""
-delay_trigger: Dict 记录哪些 reactions会引发一个delay反应：keys：reaction的idx, values: delay_trigger_affect!

-delay_complete: Dict keys:第 i 个delay channel values 完成后会引起一个 (post-delay) state-update, values: stoichiometric vector 长度是反应物长度

-delay_interrupt: Dict reactions 能够对 delay channel中造成影响的 keys : reaction idx ->  returns a Function :  how the molecules in a channel or multiple channels  will be consumed

- delay_trigger_set: 把所有会引发delay 的 reactions 的 idx 收集起来

- delay_interrupt_set: 把所有会引发delay channel 改变的 reactions 的 idx 收集起来
"""
struct DelayJumpSet
    delay_trigger::Dict
    delay_complete::Dict{Int,Vector{Pair}}
    delay_interrupt::Dict
    delay_trigger_set::Vector{Int}
    delay_interrupt_set::Vector{Int}
end
DelayJumpSet(delay_trigger,delay_complete,delay_interrupt) = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt, collect(keys(delay_trigger)), collect(keys(delay_interrupt)))

struct DelayJumpSystem
    jumpset::JumpSet
    delayjumpsets::DelayJumpSet
end
