# Switching in EMTSim: interpolation and chatter removal

EMTSim solves the network at fixed time steps Δt (Dommel's method, trapezoidal rule). Switching
does not respect that grid, which causes two separate problems. Both are described in Watson &
Arrillaga, *Power Systems Electromagnetic Transients Simulation*, section 9.4; this page says
what EMTSim does about them and where the code is.

Both methods are on by default. A circuit switches them off with `.options chatter=off interp=off`,
a run with `--chatter=off --interp=off`, and the editor has a tick box for each under Run settings.
Try them on `circuits/demo_chatter.net`, `demo_interpolation.net` and `demo_breaker.net`.

## Problem 1: a switching falls between time points

A switch time, a breaker's current zero or a diode's turn-off almost never lands exactly on a
time point. With a fixed grid the change can only be applied at the next point, so it happens up
to one step late, and a current that should have been zero is chopped.

```
 true switching ─────────┐
                         ▼
 time points ──┼────────┼────────┼────────┼──   fixed grid: the switch acts here ──┐
              t        t+Δt     t+2Δt    t+3Δt                                      ▼
```

### Method: interpolate back to the switching instant

1. Take a normal step from *t* to *t*+Δt.
2. Ask every switching element where in that step it should have switched, as a fraction
   α ∈ [0, 1] of the step (`eventAlpha` / `event_alpha` in the element classes):
   - a switch or fault with a scheduled time *t*ₛ: α = (*t*ₛ − *t*) / Δt;
   - a breaker opening at a current zero: the zero crossing of its current, found by linear
     interpolation between *i*(*t*) and *i*(*t*+Δt) (or where |*i*| drops below its chopping level);
   - a diode or thyristor turning on: where its voltage crosses the forward drop, with the gate
     fired (a thyristor whose gate fires inside the step: the firing instant, by bisection);
   - a diode or thyristor turning off: the zero crossing of its current.
3. Take the earliest α. Interpolate the **whole network state** (node voltages and every element's
   currents, voltages, fluxes, machine speeds, line history samples) linearly back to that instant
   *t*ₑ = *t* + αΔt.
4. Switch the element there, rebuild and refactorise *G*, and continue from *t*ₑ.

The solution now runs on a grid shifted by *t*ₑ. The next regular step is interpolated back onto
the original grid, so the output always has the same time points as a plain run. If several
elements switch at the same instant they are switched together; a free-wheeling diode that must
take over an inductor current turns on in the first step after the switching.

```
 ──┼──────────┼─────────────┼──          step t → t+Δt, event at α
  t          te           t+Δt
   ◄── interpolate back ──┘               state(te) = (1-α) state(t) + α state(t+Δt)
             │ switch, refactorise G
             └──── continue from te; the next step is interpolated back onto the grid
```

## Problem 2: chatter

Even when a switching lands exactly on a time point, the trapezoidal rule leaves a sustained
oscillation. For an inductor the rule is

  *i*(*t*) = *i*(*t*−Δt) + (Δt / 2*L*) [*v*(*t*) + *v*(*t*−Δt)]

Open a switch in series with it: the current is forced to zero, so the bracket must be zero at
every step, i.e. *v*(*t*) = −*v*(*t*−Δt). The inductor voltage flips sign every step, forever,
around the true value (chapter 8 of the *Tiny book of EMT simulation* shows it). Capacitor
currents do the same after a voltage jump. This is numerical, not physical.

### Method: critical damping adjustment (CDA)

After every switching, the next two steps are **backward-Euler half steps** of Δt/2. Backward
Euler has no memory of the previous voltage: for the inductor

  *i*(*t*) = *i*(*t*−Δt/2) + (Δt / 2*L*) *v*(*t*)

so the oscillation cannot start. With the half step length, its conductances (Δt/2*L*, 2*C*/Δt)
are exactly the trapezoidal ones, so *G* is not rebuilt; only the history currents change (they
simply drop their voltage term). After the two half steps, trapezoidal integration resumes.

```
 ──┼──────┼──────┼─────────────┼──   te: switching
  te   te+Δt/2  te+Δt        te+2Δt   two backward-Euler half steps, then trapezoidal again,
                                       interpolated back onto the grid
```

Why both together: interpolation alone still leaves the trapezoidal rule's half-step lag after
the switching (and the chatter); chatter removal alone acts one step late. With both, a switch
closing between time points acts at its instant to the accuracy of the trapezoidal rule itself
(see `demo_interpolation.net`).

## Where the code is, and how to change it

| Piece | Python | C++ |
|---|---|---|
| The time loop with both methods | `switching.py` → `run_with_events()` | `src/switching.cpp` → `runWithEvents()` |
| **Interpolation method** (state between two solution points) | `switching.py` → `interpolate(x, y, al)` | `src/switching.cpp` → `interpolateState()` |
| **Chatter removal method** (steps taken after a switching) | `switching.py` → `steps_after_switching()` | `src/switching.cpp` → `stepsAfterSwitching()` |
| Where each element switches within a step | `elements.py` → `event_alpha()` of `Switch`, `Diode` | `src/elements.cpp` → `eventAlpha()` |
| Backward-Euler history currents | `elements.py`, where `config.euler` is read | `src/elements.cpp`, where `config::euler` is read |
| Saving / loading / one step / earliest event | `emtcls.py` → `save`, `load`, `solve_step`, `detect` | `src/network.cpp` → `save`, `load`, `solveStep`, `detect` |
| Plain loop (both off) | `main.py` | `src/main.cpp` |

To try another chatter-removal scheme, change `steps_after_switching()`: it returns the steps to
take after a switching as `(fraction of Δt, backward Euler?)`. `[(0.5, True), (0.5, True)]` is
CDA; `[(1.0, True)]` is one full backward-Euler step (more damping); `[]` is none. Fractions need
not add up to 1: the loop interpolates back onto the output grid whatever the steps are. A
backward-Euler step of a length other than Δt/2 would need its own conductances, i.e. a change to
the elements' `stamp()` as well.

To try another interpolation, change `interpolate()` / `interpolateState()`. It receives the two
state vectors and the fraction; anything that returns a state vector works (the loop never looks
inside the vector). The state vector is the node voltages followed by every element's state in
element order (`save_state` / `saveState`).

To change where an element switches within a step (for example a quadratic fit of the current
zero instead of a linear one), change that element's `event_alpha()` / `eventAlpha()`.

Both ports implement the same algorithm and are compared on every example circuit by
`compare.py`; after a change, run it, the C++ tests (`build.bat`) and the three demos.

## Notes and limits

- The output grid is unchanged by either method; extra solution points are internal.
- The two solvers agree to round-off on all examples. Circuits with diodes or thyristors agree to
  about 10⁻⁴ of the peak, because the switching instant is the ratio of two small numbers.
- If valves keep switching each other at the same instant (more than 20 times with no progress in
  time), the step is solved the plain way, letting `check()` settle them, and the loop goes on.
- The book's one-step spike when a GTO opens between time points (its figure 9.5) does not occur
  in EMTSim even with interpolation off: a step is re-solved until the diodes agree with it.
- Chatter started by a step change of a source (not a switching) is not treated; a smaller time
  step or a snubber helps there.
