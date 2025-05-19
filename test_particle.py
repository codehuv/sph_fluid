import unittest
from math import sqrt
from particle_ import Particle # Assuming particle_.py is in the same directory or PYTHONPATH
from config import Config

# Get config values for tests
(
    N_cfg, SIM_W_cfg, BOTTOM_cfg, DAM_cfg, DAM_BREAK_cfg, G_cfg, SPACING_cfg, K_cfg, K_NEAR_cfg,
    REST_DENSITY_cfg, R_cfg, SIGMA_cfg, MAX_VEL_cfg, WALL_DAMP_cfg, VEL_DAMP_cfg, GRID_CELL_SIZE_cfg
) = Config().return_config()

class TestParticle(unittest.TestCase):

    def test_particle_initialization(self):
        p = Particle(1.0, 2.0)
        self.assertEqual(p.x_pos, 1.0)
        self.assertEqual(p.y_pos, 2.0)
        self.assertEqual(p.previous_x_pos, 1.0)
        self.assertEqual(p.previous_y_pos, 2.0)
        self.assertEqual(p.visual_x_pos, 1.0)
        self.assertEqual(p.visual_y_pos, 2.0)
        self.assertEqual(p.rho, 0.0)
        self.assertEqual(p.rho_near, 0.0)
        self.assertEqual(p.press, 0.0)
        self.assertEqual(p.press_near, 0.0)
        self.assertEqual(p.neighbors, [])
        self.assertEqual(p.x_vel, 0.0)
        self.assertEqual(p.y_vel, 0.0)
        self.assertEqual(p.x_force, 0.0)
        self.assertAlmostEqual(p.y_force, -G_cfg)

    def test_update_state_verlet_integration(self):
        p = Particle(0.0, 0.0)
        p.x_vel = 0.1
        p.y_vel = 0.2
        
        # Simulate external force calculation for this step (these are a(t))
        p.x_force = 0.05 
        p.y_force = -G_cfg + 0.01 
        
        dt = 1.0
        
        # Store initial values for comparison and calculation
        initial_x_pos = p.x_pos
        initial_y_pos = p.y_pos
        initial_x_vel = p.x_vel
        initial_y_vel = p.y_vel
        current_x_force = p.x_force # This is a(t)
        current_y_force = p.y_force # This is a(t)

        p.update_state(dam=False, dt=dt)

        # Velocity Verlet steps:
        # v(t + dt/2) = v(t) + a(t) * dt/2
        half_dt_x_vel = initial_x_vel + 0.5 * dt * current_x_force
        half_dt_y_vel = initial_y_vel + 0.5 * dt * current_y_force
        
        # x(t + dt) = x(t) + v(t + dt/2) * dt
        expected_x_pos = initial_x_pos + half_dt_x_vel * dt
        expected_y_pos = initial_y_pos + half_dt_y_vel * dt
        
        # Assume a(t+dt) is same as a(t) for this specific form of Verlet, or forces are updated externally
        # For this function, it uses current_x_force (a(t)) for the second half too.
        # v(t + dt) = v(t + dt/2) + a(t) * dt/2  (if a(t) is used for a(t+dt))
        expected_x_vel = half_dt_x_vel + 0.5 * dt * current_x_force
        expected_y_vel = half_dt_y_vel + 0.5 * dt * current_y_force
        
        self.assertAlmostEqual(p.x_pos, expected_x_pos)
        self.assertAlmostEqual(p.y_pos, expected_y_pos)
        # Velocity capping might affect these if MAX_VEL is small
        # For this test, assume MAX_VEL is large enough not to interfere
        vel_mag = sqrt(expected_x_vel**2 + expected_y_vel**2)
        if vel_mag > MAX_VEL_cfg:
            ratio = MAX_VEL_cfg / vel_mag
            expected_x_vel *= ratio
            expected_y_vel *= ratio
            
        self.assertAlmostEqual(p.x_vel, expected_x_vel)
        self.assertAlmostEqual(p.y_vel, expected_y_vel)
        
        self.assertEqual(p.previous_x_pos, initial_x_pos)
        self.assertEqual(p.previous_y_pos, initial_y_pos)
        
        # Visual position should reflect actual position unless bounded by walls
        self.assertEqual(p.visual_x_pos, p.x_pos) 
        self.assertEqual(p.visual_y_pos, p.y_pos)
        
        # Check forces are reset (assuming no wall interaction for this specific test focus)
        self.assertEqual(p.x_force, 0.0) 
        self.assertAlmostEqual(p.y_force, -G_cfg)

    def test_update_state_velocity_capping(self):
        p = Particle(0.0, 0.0)
        # Set velocity components such that magnitude is > MAX_VEL_cfg
        p.x_vel = MAX_VEL_cfg * 0.8  
        p.y_vel = MAX_VEL_cfg * 0.8 
        # Total initial velocity sqrt((0.8M)^2 + (0.8M)^2) = sqrt(1.28 M^2) = M * sqrt(1.28) approx 1.13 * M
        
        # Ensure force does not interfere significantly or cause confusion
        p.x_force = 0.0
        p.y_force = -G_cfg # Standard gravity

        p.update_state(dam=False)
        
        final_velocity_magnitude = sqrt(p.x_vel**2 + p.y_vel**2)
        # Allow for small floating point inaccuracies
        self.assertLessEqual(final_velocity_magnitude, MAX_VEL_cfg + 1e-9)

    def test_update_state_wall_constraints(self):
        # Test left wall
        p_left = Particle(-SIM_W_cfg - 0.1, BOTTOM_cfg + 0.5) # Well within y bounds
        p_left.update_state(dam=False)
        # x_force should increase due to wall: x_force -= 0.3 * (x_pos - -SIM_W) * WALL_DAMP
        # x_pos - -SIM_W is negative, so -= (negative) means x_force increases (pushed right)
        self.assertGreater(p_left.x_force, 0.0 + 1e-9) # Initial x_force is 0, wall adds positive
        self.assertEqual(p_left.visual_x_pos, -SIM_W_cfg)

        # Test right wall (no dam)
        p_right = Particle(SIM_W_cfg + 0.1, BOTTOM_cfg + 0.5)
        p_right.update_state(dam=False)
        # x_force should decrease due to wall: x_force -= 0.3 * (x_pos - SIM_W) * WALL_DAMP
        # x_pos - SIM_W is positive, so -= (positive) means x_force decreases (pushed left)
        self.assertLess(p_right.x_force, 0.0 - 1e-9) # Initial x_force is 0, wall adds negative
        self.assertEqual(p_right.visual_x_pos, SIM_W_cfg)

        # Test dam wall
        p_dam = Particle(DAM_cfg + 0.1, BOTTOM_cfg + 0.5) # Assuming DAM_cfg < SIM_W_cfg
        p_dam.update_state(dam=True) # Dam is active
        # x_force should decrease: x_force -= (x_pos - DAM) * WALL_DAMP
        self.assertLess(p_dam.x_force, 0.0 - 1e-9) 
        # No visual clamping for dam in current code, visual_x_pos will be actual x_pos

        # Test bottom wall
        p_bottom = Particle(0.0, BOTTOM_cfg - 0.1) # Well within x bounds
        initial_y_force_before_wall_effect = -G_cfg # Force before wall constraint is applied in update_state
        
        p_bottom.update_state(dam=False)
        # y_force -= 0.7 * (y_pos - SIM_W_cfg) * WALL_DAMP
        # (y_pos - SIM_W_cfg) = (BOTTOM_cfg - 0.1 - SIM_W_cfg) which is very negative.
        # So, y_force -= (very_negative_val), meaning y_force increases (pushed up).
        # The final p.y_force is -G_cfg (from reset) + wall_push_y
        # So, if wall_push_y is positive, p.y_force should be > -G_cfg
        self.assertGreater(p_bottom.y_force, -G_cfg - (0.7 * ( (BOTTOM_cfg - 0.1) - SIM_W_cfg) * WALL_DAMP_cfg) -1e-9 )
        self.assertEqual(p_bottom.visual_y_pos, BOTTOM_cfg)


    def test_calculate_pressure(self):
        p = Particle(0.0, 0.0)
        p.rho = REST_DENSITY_cfg + 1.0 # Example rho
        p.rho_near = 2.0 # Example rho_near
        
        p.calculate_pressure()
        
        expected_press = K_cfg * (p.rho - REST_DENSITY_cfg)
        expected_press_near = K_NEAR_cfg * p.rho_near
        
        self.assertAlmostEqual(p.press, expected_press)
        self.assertAlmostEqual(p.press_near, expected_press_near)

if __name__ == '__main__':
    unittest.main()
