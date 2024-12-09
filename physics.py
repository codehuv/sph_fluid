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


def calculate_density(particles: list[Particle]) -> None:
    """
    입자의 밀도를 계산합니다.
        밀도는 이웃 입자 간의 상대 거리를 합산하여 계산됩니다.
        입자가 서로 충돌하여 불안정성을 야기하는 것을 방지하기 위해
        밀도와 근접 밀도를 구분하여 계산합니다.

    Args:
        particles (list[Particle]): 입자 리스트
    """
    for i, particle_1 in enumerate(particles):
        density = 0.0
        density_near = 0.0
        # 주변 이웃간 밀도 계산
        for particle_2 in particles[i + 1 :]:
            distance = sqrt(
                (particle_1.x_pos - particle_2.x_pos) ** 2
                + (particle_1.y_pos - particle_2.y_pos) ** 2
            )
            if distance < R:
                # 기본 거리 1로 설정
                normal_distance = 1 - distance / R
                density += normal_distance**2
                density_near += normal_distance**3
                particle_2.rho += normal_distance**2
                particle_2.rho_near += normal_distance**3
                particle_1.neighbors.append(particle_2)
        particle_1.rho += density
        particle_1.rho_near += density_near


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