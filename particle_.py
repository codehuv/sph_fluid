from math import sqrt
from config import Config


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


class Particle:
    """
    속성:
    x_pos: 입자의 x 위치
    y_pos: 입자의 y 위치
    previous_x_pos: 이전 프레임에서 입자의 x 위치
    previous_y_pos: 이전 프레임에서 입자의 y 위치
    visual_x_pos: 화면에 표시되는 입자의 x 위치
    visual_y_pos: 화면에 표시되는 입자의 y 위치
    rho: 입자의 밀도
    rho_near: 입자의 근접 밀도, 입자 간 충돌 방지에 사용됨
    press: 입자의 압력
    press_near: 입자의 근접 압력, 입자 간 충돌 방지에 사용됨
    neighbors: 입자의 이웃 입자 목록
    x_vel: 입자의 x 속도
    y_vel: 입자의 y 속도
    x_force: 입자에 가해지는 x 방향 힘
    y_force: 입자에 가해지는 y 방향 힘
    """

    def __init__(self, x_pos: float, y_pos: float):
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.previous_x_pos = x_pos
        self.previous_y_pos = y_pos
        self.visual_x_pos = x_pos
        self.visual_y_pos = y_pos
        self.rho = 0.0
        self.rho_near = 0.0
        self.press = 0.0
        self.press_near = 0.0
        self.neighbors = []
        self.x_vel = 0.0
        self.y_vel = 0.0
        self.x_force = 0.0
        self.y_force = -G

    def update_state(self, dam: bool, dt: float = 1.0):
        """
        Updates the particle's state using the Velocity Verlet integration method.

        Args:
            dam (bool): Indicates whether the dam is present.
            dt (float, optional): The time step. Defaults to 1.0.
        """

        # 이전 위치 보존
        self.previous_x_pos = self.x_pos
        self.previous_y_pos = self.y_pos

        # 1. 이전 속도와 현재 힘을 이용하여 속도의 절반 단계를 계산 (half-step velocity)
        half_x_vel = self.x_vel + 0.5 * dt * self.x_force
        half_y_vel = self.y_vel + 0.5 * dt * self.y_force

        # 2. 절반 단계의 속도를 사용하여 위치 업데이트
        self.x_pos += half_x_vel * dt
        self.y_pos += half_y_vel * dt

        # 3. 새로운 위치에서의 힘 계산은 바깥(main.py)에서 해 줄 예정

        # 4. 새로운 위치에서 계산된 힘을 사용하여 속도 업데이트
        self.x_vel = half_x_vel + 0.5 * dt * self.x_force
        self.y_vel = half_y_vel + 0.5 * dt * self.y_force

        # 화면에 표시되는 시각적 위치 설정
        self.visual_x_pos = self.x_pos
        self.visual_y_pos = self.y_pos
        
        # force 초기화
        (self.x_force, self.y_force) = (0.0, -G)

        # 속도 계산 (Verlet에서는 덜 중요하지만, 필요에 따라 계산)
        velocity = sqrt(self.x_vel**2 + self.y_vel**2)

        # 속도가 너무 높으면 감소시킴
        if velocity > MAX_VEL:
            reduction_ratio = MAX_VEL / velocity
            self.x_vel *= reduction_ratio
            self.y_vel *= reduction_ratio

        # 벽 제약 조건
        if self.x_pos < -SIM_W:
            self.x_force -= 0.3 * (self.x_pos - -SIM_W) * WALL_DAMP
            self.visual_x_pos = -SIM_W
        if dam is True and self.x_pos > DAM:
            self.x_force -= (self.x_pos - DAM) * WALL_DAMP
        if self.x_pos > SIM_W:
            self.x_force -= 0.3 * (self.x_pos - SIM_W) * WALL_DAMP
            self.visual_x_pos = SIM_W
        if self.y_pos < BOTTOM:
            self.y_force -= 0.7 * (self.y_pos - SIM_W) * WALL_DAMP
            self.visual_y_pos = BOTTOM

        # 밀도 초기화
        self.rho = 0.0
        self.rho_near = 0.0

        # 이웃 입자 목록 초기화
        self.neighbors = []

    def calculate_pressure(self):
        """
        입자의 압력을 계산
        """
        self.press = K * (self.rho - REST_DENSITY)
        self.press_near = K_NEAR * self.rho_near