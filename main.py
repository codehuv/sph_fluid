import numpy as np
from config import Config
from particle_ import Particle
import pygame
from physics import (
    start,
    calculate_density,
    create_pressure,
    calculate_viscosity,
    create_grid
)

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

def update(particles: list[Particle], dam: bool) -> list[Particle]:
    """
    Calculates one step of the simulation.
    """
    # 1. 힘 초기화, 중력가속도 적용, 벽과의 상호작용
    for particle in particles:
        particle.force = np.array([0.0, -G], dtype=float)
        #particle.update_state(dam) # 여기에서 불리면 안됨

    # 2. 밀도 계산
    grid = create_grid(particles, GRID_CELL_SIZE)
    calculate_density(particles, grid, GRID_CELL_SIZE)

    # 3. 압력 계산
    for particle in particles:
        particle.calculate_pressure()

    # 4. 압력 힘 적용
    create_pressure(particles)

    # 5. 점성 힘 적용
    calculate_viscosity(particles)

    # 6. 업데이트된 힘을 바탕으로 update_state 호출
    for particle in particles:
        particle.update_state(dam)

    return particles

# Pygame 설정
pygame.init()
screen_width = int(2.75 * SIM_W * 100)
screen_height = int(1 * SIM_W * 100)
screen = pygame.display.set_mode((screen_width, screen_height))
pygame.display.set_caption("2D SPH particle interaction simulation")
particle_radius = int(SPACING * 50)  # 필요에 따라 입자 크기 조정

simulation_state = start(-SIM_W, SIM_W, BOTTOM+1, 0.03, N)

frame = 0
dam_built = False
running = True

# 시뮬레이션 좌표를 화면 좌표로 변환하는 함수
def sim_to_screen(x, y):
    screen_x = int((x + SIM_W) * 100)
    screen_y = int((SIM_W - y) * 100)  # Pygame 좌표계에 맞게 y축 반전
    return screen_x, screen_y

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    simulation_state = update(simulation_state, dam_built)

    # 화면 지우기
    screen.fill((0, 0, 0))

    # 입자 그리기
    for particle in simulation_state:
        screen_x, screen_y = sim_to_screen(particle.visual_x_pos, particle.visual_y_pos)
        pygame.draw.circle(screen, (0, 0, 255), (screen_x, screen_y), particle_radius)

    # 화면 업데이트
    pygame.display.flip()

    frame += 1
    pygame.time.delay(5)  # 애니메이션 속도에 맞게 지연 시간 조정

pygame.quit()