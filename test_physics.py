import unittest
from math import sqrt
from particle_ import Particle # Assuming particle_.py is in the same directory or PYTHONPATH
from physics import start, calculate_density, create_pressure, calculate_viscosity, create_grid
from config import Config

# Get config values
(
    N_cfg, SIM_W_cfg, BOTTOM_cfg, DAM_cfg, DAM_BREAK_cfg, G_cfg, SPACING_cfg, K_cfg, K_NEAR_cfg,
    REST_DENSITY_cfg, R_cfg, SIGMA_cfg, MAX_VEL_cfg, WALL_DAMP_cfg, VEL_DAMP_cfg, GRID_CELL_SIZE_cfg
) = Config().return_config()

class TestPhysicsFunctions(unittest.TestCase):

    def test_start_particle_creation(self):
        # Test with specific parameters where behavior of `xmax-1` is predictable
        xmin, ymin = 0.0, 0.0
        space = 0.1
        
        # Case 1: Fill less than one row
        count1 = 2
        xmax1 = 0.3 # xmax-1 = -0.7. x_pos=0.0, 0.1, 0.2. Next x_pos=0.3. 0.3 > -0.7.
                    # The condition `if x_pos > xmax-1` means wrapping happens if new x_pos > xmax-1
        particles1 = start(xmin, xmax1, ymin, space, count1)
        self.assertEqual(len(particles1), count1)
        self.assertAlmostEqual(particles1[0].x_pos, 0.0)
        self.assertAlmostEqual(particles1[0].y_pos, 0.0)
        self.assertAlmostEqual(particles1[1].x_pos, 0.1)
        self.assertAlmostEqual(particles1[1].y_pos, 0.0)

        # Case 2: Force wrapping
        count2 = 3
        xmax2 = 0.15 # xmax-1 = -0.85.
                     # p0: (0,0). x_pos=0.1. 0.1 > -0.85 (No wrap)
                     # p1: (0.1,0). x_pos=0.2. 0.2 > -0.85 (No wrap based on condition as is)
                     # This `xmax-1` logic is peculiar. Let's assume `xmax` is a boundary for placement.
                     # If xmax = 0.15, particles at 0.0, 0.1. next x_pos = 0.2.
                     # The actual condition `if x_pos > xmax-1` means:
                     # For xmax2=0.15, xmax2-1 = -0.85
                     # p0: x=0.0. next x_pos=0.1. (0.1 > -0.85) is true, but this doesn't mean wrap.
                     # It means x_pos is incremented. It wraps *after* adding particle *if* condition is met.
                     # x_pos = xmin (0.0) -> p0(0.0, 0.0). x_pos becomes 0.1. Is 0.1 > -0.85? Yes.
                     # x_pos = 0.1 -> p1(0.1, 0.0). x_pos becomes 0.2. Is 0.2 > -0.85? Yes.
                     # The logic seems to be: place particle, then update x_pos. If new x_pos is too large, wrap.
                     # If xmax is the limit of the particle positions:
                     # Particle at 0.0. next_x = 0.1.
                     # Particle at 0.1. next_x = 0.2. If 0.2 is considered "out of bounds for the row", then wrap.
                     # The current condition `x_pos > xmax-1` for xmax=0.15 means `x_pos > -0.85`.
                     # This will always be true for positive x_pos.
                     # This implies the number of particles per row is ((xmax-1) - xmin) / space.
                     # This needs re-evaluation or testing with the simulation's actual parameters.
                     # Let's use xmax that ensures some particles per row before wrapping.
        xmax3 = 0.2 # xmax-1 = -0.8. (Still problematic for small positive x_pos)
                    # If xmax is effectively SIM_W (e.g., 3), then xmax-1=2.
                    # xmin=0, xmax_boundary=0.2, space=0.1, count=3
                    # p0: (0.0, 0.0). x_pos becomes 0.1.
                    # p1: (0.1, 0.0). x_pos becomes 0.2.
                    # if x_pos (0.2) > xmax_boundary-1 (-0.8), it wraps.
                    # This means positions: (0,0), (0.1,0), then (0, 0.1)
        particles3 = start(xmin=0.0, xmax=0.2, ymin=0.0, space=0.1, count=3)
        self.assertEqual(len(particles3), 3)
        self.assertAlmostEqual(particles3[0].x_pos, 0.0)
        self.assertAlmostEqual(particles3[0].y_pos, 0.0)
        self.assertAlmostEqual(particles3[1].x_pos, 0.1)
        self.assertAlmostEqual(particles3[1].y_pos, 0.0)
        self.assertAlmostEqual(particles3[2].x_pos, 0.0) # Wrapped
        self.assertAlmostEqual(particles3[2].y_pos, 0.1)


    def test_create_grid(self):
        # Config: SIM_W_cfg = 3, R_cfg = 0.08*1.25 = 0.1, GRID_CELL_SIZE_cfg = 0.1*1.5 = 0.15
        p1 = Particle(0.0, 0.0) 
        p2 = Particle(0.1, 0.0) # Should be in the same cell as p1 or adjacent
        p3 = Particle(0.2, 0.0) # Should be in a different cell from p1 if GRID_CELL_SIZE is small enough
        particles = [p1, p2, p3]
        
        grid = create_grid(particles, GRID_CELL_SIZE_cfg)

        # Cell index: int((pos + SIM_W_cfg) / GRID_CELL_SIZE_cfg)
        # p1 (0,0): cell_x = int((0+3)/0.15) = int(20) = 20. cell_y = int((0+3)/0.15) = 20. Key (20,20)
        # p2 (0.1,0): cell_x = int((0.1+3)/0.15) = int(3.1/0.15) = int(20.66) = 20. Key (20,20)
        # p3 (0.2,0): cell_x = int((0.2+3)/0.15) = int(3.2/0.15) = int(21.33) = 21. Key (21,20)

        cell_p1_key = (int((p1.x_pos + SIM_W_cfg) / GRID_CELL_SIZE_cfg), int((p1.y_pos + SIM_W_cfg) / GRID_CELL_SIZE_cfg))
        cell_p2_key = (int((p2.x_pos + SIM_W_cfg) / GRID_CELL_SIZE_cfg), int((p2.y_pos + SIM_W_cfg) / GRID_CELL_SIZE_cfg))
        cell_p3_key = (int((p3.x_pos + SIM_W_cfg) / GRID_CELL_SIZE_cfg), int((p3.y_pos + SIM_W_cfg) / GRID_CELL_SIZE_cfg))
        
        self.assertEqual(cell_p1_key, (20,20))
        self.assertEqual(cell_p2_key, (20,20))
        self.assertEqual(cell_p3_key, (21,20))

        self.assertIn(p1, grid[cell_p1_key])
        self.assertIn(p2, grid[cell_p2_key]) # Same cell as p1
        self.assertIn(p3, grid[cell_p3_key])
        
        self.assertEqual(len(grid[cell_p1_key]), 2) # p1 and p2
        self.assertEqual(len(grid[cell_p3_key]), 1) # p3

    def test_calculate_density(self):
        # Config: R_cfg = 0.1, GRID_CELL_SIZE_cfg = 0.15
        p1 = Particle(0.0, 0.0)
        p2 = Particle(R_cfg * 0.5, 0.0) # distance = 0.05 (neighbor)
        p3 = Particle(R_cfg * 2.0, 0.0) # distance = 0.2 (not a neighbor)
        particles = [p1, p2, p3]

        grid = create_grid(particles, GRID_CELL_SIZE_cfg)
        calculate_density(particles, grid, GRID_CELL_SIZE_cfg)

        # For p1 and p2:
        dist_p1_p2 = R_cfg * 0.5
        norm_dist_factor = 1 - dist_p1_p2 / R_cfg # 1 - 0.5 = 0.5

        # p1's density from p2 (and vice-versa due to symmetry)
        expected_rho_contrib = norm_dist_factor**2      # (0.5)^2 = 0.25
        expected_rho_near_contrib = norm_dist_factor**3 # (0.5)^3 = 0.125
        
        self.assertAlmostEqual(p1.rho, expected_rho_contrib)
        self.assertAlmostEqual(p1.rho_near, expected_rho_near_contrib)
        self.assertIn(p2, p1.neighbors)
        self.assertNotIn(p3, p1.neighbors)
        self.assertEqual(len(p1.neighbors), 1)

        self.assertAlmostEqual(p2.rho, expected_rho_contrib) # Symmetric calculation for p2 from p1
        self.assertAlmostEqual(p2.rho_near, expected_rho_near_contrib)
        self.assertIn(p1, p2.neighbors)
        self.assertNotIn(p3, p2.neighbors)
        self.assertEqual(len(p2.neighbors), 1)

        # p3's density (should be 0 as p1,p2 are > R away)
        self.assertAlmostEqual(p3.rho, 0.0)
        self.assertAlmostEqual(p3.rho_near, 0.0)
        self.assertEqual(len(p3.neighbors), 0)

    def test_create_pressure(self):
        # Config: R_cfg = 0.1, K_cfg, K_NEAR_cfg, REST_DENSITY_cfg = 3.0
        p1 = Particle(0.0, 0.0)
        p2 = Particle(R_cfg * 0.5, 0.0) # dist = 0.05
        particles = [p1, p2]

        # Setup: p1 and p2 are neighbors
        p1.neighbors = [p2]
        p2.neighbors = [p1]

        # Assign densities and calculate initial pressures to create a repulsive force
        # (rho > REST_DENSITY_cfg)
        dist_p1_p2 = R_cfg * 0.5
        norm_dist_q = 1 - dist_p1_p2 / R_cfg # 0.5

        # To simplify, let's assume they have some pressure values
        # Say rho makes (rho - REST_DENSITY) = 1.0 for `press`
        # And rho_near makes rho_near = 1.0 for `press_near`
        p1.press = K_cfg * 1.0
        p1.press_near = K_NEAR_cfg * 1.0
        p2.press = K_cfg * 1.0
        p2.press_near = K_NEAR_cfg * 1.0
        
        # Store initial forces (these are typically -G in y, 0 in x)
        p1_initial_fx, p1_initial_fy = p1.x_force, p1.y_force
        p2_initial_fx, p2_initial_fy = p2.x_force, p2.y_force

        create_pressure(particles)

        # Expected behavior: p1 pushed left (negative_dx_force), p2 pushed right (positive_dx_force)
        self.assertLess(p1.x_force, p1_initial_fx)
        self.assertGreater(p2.x_force, p2_initial_fx)
        
        # Y forces should ideally be unchanged by this horizontal interaction
        self.assertAlmostEqual(p1.y_force, p1_initial_fy)
        self.assertAlmostEqual(p2.y_force, p2_initial_fy)

        # Check Newton's third law for the change in force
        delta_p1_fx = p1.x_force - p1_initial_fx
        delta_p2_fx = p2.x_force - p2_initial_fx
        self.assertAlmostEqual(delta_p1_fx, -delta_p2_fx)

    def test_calculate_viscosity(self):
        # Config: R_cfg = 0.1, SIGMA_cfg = 2.5
        p1 = Particle(0.0, 0.0)
        p2 = Particle(R_cfg * 0.5, 0.0) # dist = 0.05
        particles = [p1, p2] # Function iterates this list

        # p1 and p2 are neighbors
        p1.neighbors = [p2]
        p2.neighbors = [p1] # For full interaction if function logic depends on it

        # Case 1: Particles approaching each other
        p1.x_vel = 0.2; p1.y_vel = 0.0
        p2.x_vel = -0.2; p2.y_vel = 0.0
        
        initial_p1_x_vel, initial_p2_x_vel = p1.x_vel, p2.x_vel
        
        calculate_viscosity(particles) # Call with both particles

        # Calculation for interaction p1-p2:
        # particle_to_neighbor = [0.05, 0]. distance = 0.05. normal_p_to_n = [1, 0]
        # relative_distance (q_rel) = 0.05 / 0.1 = 0.5
        # velocity_difference_proj = ( (p1.x_vel - p2.x_vel)*normal_p_to_n[0] + ... )
        #                            = ( (0.2 - (-0.2)) * 1 ) = 0.4
        # This is > 0, so viscosity applies.
        # visc_force_magnitude_comp = (1 - q_rel) * SIGMA_cfg * velocity_difference_proj
        #                           = (1 - 0.5) * 2.5 * 0.4 = 0.5 * 2.5 * 0.4 = 0.5 * 1.0 = 0.5
        # viscosity_force_vector = [visc_force_magnitude_comp * normal_p_to_n[0], ...] = [0.5, 0]
        
        # p1.x_vel -= viscosity_force_vector[0] * 0.5
        # p2.x_vel += viscosity_force_vector[0] * 0.5
        # This happens for p1->p2 interaction.
        # If p2->p1 interaction is also processed (it will be because p2 is in `particles`):
        #   normal_p_to_n becomes normal_p2_to_p1 = [-1, 0]
        #   velocity_difference_proj = ( (p2.x_vel - p1.x_vel)*normal_p2_to_p1[0] )
        #                            = ( (-0.2 - 0.2) * -1 ) = (-0.4 * -1) = 0.4 (same magnitude)
        #   visc_force_vector for p2->p1 = [0.5 * (-1), 0] = [-0.5, 0]
        #   p2.x_vel -= (-0.5) * 0.5
        #   p1.x_vel += (-0.5) * 0.5
        # Total change for p1.x_vel: -0.5*0.5 (from p1->p2) + (-0.5)*0.5 (from p2->p1) = -0.25 - 0.25 = -0.5
        # Total change for p2.x_vel: +0.5*0.5 (from p1->p2) + (+0.5)*0.5 (from p2->p1) = +0.25 + 0.25 = +0.5
        
        expected_total_change_p1 = -0.5 
        expected_total_change_p2 = 0.5
        
        self.assertAlmostEqual(p1.x_vel, initial_p1_x_vel + expected_total_change_p1)
        self.assertAlmostEqual(p2.x_vel, initial_p2_x_vel + expected_total_change_p2)
        self.assertAlmostEqual(p1.y_vel, 0.0) # No y-component change
        self.assertAlmostEqual(p2.y_vel, 0.0)

        # Case 2: Particles moving apart (velocity_difference_proj < 0)
        p1.x_vel = -0.1; p2.x_vel = 0.1
        initial_p1_x_vel, initial_p2_x_vel = p1.x_vel, p2.x_vel
        
        calculate_viscosity(particles)
        # velocity_difference_proj for p1->p2 = ((-0.1 - 0.1)*1) = -0.2. Should not apply.
        self.assertAlmostEqual(p1.x_vel, initial_p1_x_vel)
        self.assertAlmostEqual(p2.x_vel, initial_p2_x_vel)


if __name__ == '__main__':
    unittest.main()
