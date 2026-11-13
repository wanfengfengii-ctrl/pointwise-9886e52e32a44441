"""Recorded scheduler semantics, not an inference implementation."""
from functools import lru_cache


class Design:
    def __init__(self, protocol):
        self.p = protocol
        self.n = protocol["n"]
        self.block = protocol["block"]
        self.quota = protocol["on_per_block"]
        self.hot = tuple(protocol["hot"])
        self.hot_on = protocol["hot_on"]
        self.max_run = protocol["max_run"]
        if len(self.hot) != self.n or self.n % self.block:
            raise ValueError("PROTOCOL_SHAPE")

    def step(self, state, bit):
        i, last, run, used, hot_used = state
        if i >= self.n or bit not in (0, 1):
            return None
        run = run + 1 if bit == last else 1
        used += bit
        hot_used += bit * self.hot[i]
        if run > self.max_run or used > self.quota or hot_used > self.hot_on:
            return None
        if (i + 1) % self.block == 0:
            if used != self.quota:
                return None
            used = 0
        return (i + 1, bit, run, used, hot_used)

    @lru_cache(None)
    def possible(self, state):
        if state[0] == self.n:
            return state[3] == 0 and state[4] == self.hot_on
        for bit in (0, 1):
            nxt = self.step(state, bit)
            if nxt is not None and self.possible(nxt):
                return True
        return False

    @property
    def initial(self):
        return (0, -1, 0, 0, 0)

    def allowed(self, state):
        return tuple(bit for bit in (0, 1)
                     if (nxt := self.step(state, bit)) is not None
                     and self.possible(nxt))

    def weights(self, state, previous_outcomes):
        choices = self.allowed(state)
        if len(choices) == 1:
            return {choices[0]: 3}
        if len(choices) != 2:
            raise ValueError("NO_FEASIBLE_CONTINUATION")
        one = 2 if sum(previous_outcomes[-3:]) >= self.p["threshold"] else 1
        return {0: 3 - one, 1: one}

    def path_mass(self, bits, baseline, tau, gamma):
        """Integer numerator over 3**n; only for a single specified path."""
        if len(bits) != self.n or len(baseline) != self.n:
            raise ValueError("PATH_LENGTH")
        state, mass, lag = self.initial, 1, 0
        history = list(self.p["pre_y"])
        for i, bit in enumerate(bits):
            weights = self.weights(state, history)
            if bit not in weights:
                return 0
            mass *= weights[bit]
            history.append(baseline[i] + tau * bit + gamma * lag)
            history = history[-3:]
            lag = bit
            state = self.step(state, bit)
        return mass

    def tiny_paths(self):
        """Exhaustive reference for small regressions; never for the archive."""
        if self.n > 10:
            raise ValueError("TINY_ONLY")
        def visit(state, path):
            if state[0] == self.n:
                yield tuple(path)
                return
            for bit in self.allowed(state):
                yield from visit(self.step(state, bit), path + [bit])
        yield from visit(self.initial, [])
