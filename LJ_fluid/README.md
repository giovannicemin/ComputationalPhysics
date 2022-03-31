# The Lennard-Jones Fluid
Even though this project is quite simple, I really enjoyed 
the first time I implemented the simulation. 
Dear reader, I hope you will have fun too.

**Description**: The program set up a Lennard-Jones fluid simulation
using molecular dynamics techniques.
After a phase of equilibration I procede to do some measurements on 
the fluid, namely measuring mean potential and kinetic energy and 
calculating the radial distribution function.

----

## How it works

The fluid is completely classical, so the particle's trajectories are 
calculated by Newton's second law. The algorithm uto integrate such 
equation is Velocity-Verlet:
<pre><code>
x<sub>n+1</sub> = x<sub>n</sub> + v<sub>n</sub> Δt + 1/2 a<sub>n</sub> Δt<sup>2</sup> 

</code></pre>

### Initial state

### Equilibration

### Measurements
----

## Credits and references

1. Projects that inspired you
2. Related projects
3. Books, papers, talks, or other sources that have meaningful impact or influence on this project
