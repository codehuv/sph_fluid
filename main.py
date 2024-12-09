import numpy as np
from config import Config
from particle_ import Particle
import pygame
from physics import (
    start,
    calculate_density,
    create_pressure,
    calculate_viscosity,
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
) = Config().return_config()

def update(particles: list[Particle], dam: bool) -> list[Particle]:
    """
    시뮬레이션의 한 단계를 계산합니다.
    """
    # 입자 상태 업데이트 (힘 적용, 값 초기화 등)
    for particle in particles:
        particle.update_state(dam)

    # 밀도 계산
    calculate_density(particles)

    # 압력 계산
    for particle in particles:
        particle.calculate_pressure()

    # 압력 힘 적용
    create_pressure(particles)

    # 점성 힘 적용
    calculate_viscosity(particles)

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