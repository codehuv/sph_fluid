from math import sqrt

from config import Config
from particle_ import Particle


(
    N,
    SIM_W,
    BOTTOM,
    DAM,
    DAM_BREAK,
    G,
    SPACING,
    K,
    K_NEAR,
    REST_DENSITY,
    R,
    SIGMA,
    MAX_VEL,
    WALL_DAMP,
    VEL_DAMP,
    GRID_CELL_SIZE
) = Config().return_config()


def start(
    xmin: float, xmax: float, ymin: float, space: float, count: int
) -> list[Particle]:
    """
    xmin, xmax, ymin 범위 내에 입자 사각형을 생성합니다.
    (xmin, ymin) 위치에서 입자를 생성하고,
    입자 수가 count에 도달할 때까지 입자를 추가합니다.
    입자는 위치 [x, y]로 표현됩니다.

    Args:
        xmin (float): 사각형의 x 최소 경계
        xmax (float): 사각형의 x 최대 경계
        ymin (float): 사각형의 y 최소 경계
        space (float): 입자 간 간격
        count (int): 입자 수

    Returns:
        list: Particle 객체 리스트
    """
    result = []
    x_pos, y_pos = xmin, ymin
    for _ in range(count):
        result.append(Particle(x_pos, y_pos))
        x_pos += space
        if x_pos > xmax-1:
            x_pos = xmin
            y_pos += space
    return result


def calculate_density(particles: list[Particle], grid: dict, grid_cell_size: float) -> None:
    """
    Calculates the density and near-density of each particle.

    Args:
        particles (list[Particle]): The list of particles.
        grid (dict): The grid containing particles assigned to cells.
        grid_cell_size (float): The size of each grid cell.
    """
    for particle in particles:
        particle.rho = 0.0
        particle.rho_near = 0.0
        particle.neighbors = []

        cell_x = int((particle.x_pos + SIM_W) / grid_cell_size)
        cell_y = int((particle.y_pos + SIM_W) / grid_cell_size)

        # Iterate through neighboring cells (including the particle's own cell)
        for i in range(cell_x - 1, cell_x + 2):
            for j in range(cell_y - 1, cell_y + 2):
                if (i, j) in grid:
                    for neighbor in grid[(i, j)]:
                        if particle != neighbor:  # Exclude self
                            distance = sqrt(
                                (particle.x_pos - neighbor.x_pos) ** 2
                                + (particle.y_pos - neighbor.y_pos) ** 2
                            )
                            if distance < R:
                                normal_distance = 1 - distance / R
                                particle.rho += normal_distance**2
                                particle.rho_near += normal_distance**3
                                particle.neighbors.append(neighbor)


def create_pressure(particles: list[Particle]) -> None:
    """
    입자의 압력 힘을 계산합니다.
        calculate_density 함수에서 이웃 리스트와 압력이 이미 계산됨
        각 이웃 입자의 압력 힘을 합산하여 압력 힘을 계산하고,
        이웃 입자 방향으로 힘을 적용합니다.

    Args:
        particles (list[Particle]): 입자 리스트
    """
    for particle in particles:
        press_x = 0.0
        press_y = 0.0
        for neighbor in particle.neighbors:
            particle_to_neighbor = [
                neighbor.x_pos - particle.x_pos,
                neighbor.y_pos - particle.y_pos,
            ]
            distance = sqrt(particle_to_neighbor[0] ** 2 + particle_to_neighbor[1] ** 2)
            normal_distance = 1 - distance / R
            total_pressure = (
                particle.press + neighbor.press
            ) * normal_distance**2 + (
                particle.press_near + neighbor.press_near
            ) * normal_distance**3
            pressure_vector = [
                particle_to_neighbor[0] * total_pressure / distance,
                particle_to_neighbor[1] * total_pressure / distance,
            ]
            neighbor.x_force += pressure_vector[0]
            neighbor.y_force += pressure_vector[1]
            press_x += pressure_vector[0]
            press_y += pressure_vector[1]
        particle.x_force -= press_x
        particle.y_force -= press_y


def calculate_viscosity(particles: list[Particle]) -> None:
    """
    입자의 점성 힘을 계산합니다.
    힘 = (입자 간 상대 거리) * (점성 가중치) * (입자 간 속도 차이)
    속도 차이는 입자 사이의 벡터를 기반으로 계산됩니다.

    Args:
        particles (list[Particle]): 입자 리스트
    """

    for particle in particles:
        for neighbor in particle.neighbors:
            particle_to_neighbor = [
                neighbor.x_pos - particle.x_pos,
                neighbor.y_pos - particle.y_pos,
            ]
            distance = sqrt(particle_to_neighbor[0] ** 2 + particle_to_neighbor[1] ** 2)
            normal_p_to_n = [
                particle_to_neighbor[0] / distance,
                particle_to_neighbor[1] / distance,
            ]
            relative_distance = distance / R
            velocity_difference = (particle.x_vel - neighbor.x_vel) * normal_p_to_n[
                0
            ] + (particle.y_vel - neighbor.y_vel) * normal_p_to_n[1]
            if velocity_difference > 0:
                viscosity_force = [
                    (1 - relative_distance)
                    * SIGMA
                    * velocity_difference
                    * normal_p_to_n[0],
                    (1 - relative_distance)
                    * SIGMA
                    * velocity_difference
                    * normal_p_to_n[1],
                ]
                particle.x_vel -= viscosity_force[0] * 0.5
                particle.y_vel -= viscosity_force[1] * 0.5
                neighbor.x_vel += viscosity_force[0] * 0.5
                neighbor.y_vel += viscosity_force[1] * 0.5

def create_grid(particles: list[Particle], grid_cell_size: float) -> dict:
    grid = {}
    for particle in particles:
        cell_x = int((particle.x_pos + SIM_W) / grid_cell_size)
        cell_y = int((particle.y_pos + SIM_W) / grid_cell_size)
        cell_coords = (cell_x, cell_y)
        if cell_coords not in grid:
            grid[cell_coords] = []
        grid[cell_coords].append(particle)
    return grid