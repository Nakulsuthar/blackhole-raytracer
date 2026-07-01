import numpy as np
import sys
sys.path.insert(0, sys.path[0] + '/..')
from tracer.integrator import RayTracer

rt = RayTracer(rs=1.0, r0 = 50.0)

# Test 1: ray below b_crit should be captured
captured, phi = rt.trace_ray(rt.b_crit * 0.9)
assert captured, "FAIL: ray at 0.9 * b_crit should be captured"
print(f"Test 1 PASS — b = 0.9 * b_crit: captured at phi = {phi:.3f}")

# Test 2: ray above b_crit should escape
captured, phi = rt.trace_ray(rt.b_crit * 1.1)
assert not captured, "FAIL: ray at 1.1 * b_crit should escape"
print(f"Test 2 PASS — b = 1.1 * b_crit: escaped at phi = {phi:.3f}")
# Test 3: weak-field deflection
b_test = 20.0 * rt.rs   # larger b = more firmly in weak field
captured, exit_phi = rt.trace_ray(b_test)
assert not captured, "FAIL: ray at b = 20 rs should escape"

# angle a straight line would sweep at finite r0
flat_angle = np.pi - 2 * np.arcsin(b_test / rt.r0)

# deflection is the excess beyond the flat-space angle
deflection = exit_phi - flat_angle
expected = 2 * rt.rs / b_test

error = abs(deflection - expected) / expected
assert error < 0.10, f"FAIL: deflection {deflection:.6f} vs expected {expected:.6f} (error {error:.1%})"
print(f"Test 3 PASS — weak field: deflection = {deflection:.6f}, expected = {expected:.6f}, error = {error:.1%}")
print("\nAll tests passed.")