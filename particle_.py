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

    def update_state(self, dam: bool):
        # 이전의 위치 초기화
        (self.previous_x_pos, self.previous_y_pos) = (self.x_pos, self.y_pos)

        # 뉴턴의 제2 법칙과 오일러 적분을 사용하여 힘 적용 (질량 = 1, dt = 1)
        (self.x_vel, self.y_vel) = (
            self.x_vel + self.x_force,
            self.y_vel + self.y_force,
        )

        # 오일러 적분(dt = 1)을 사용하여 속도에 따라 입자 이동
        (self.x_pos, self.y_pos) = (self.x_pos + self.x_vel, self.y_pos + self.y_vel)

        # 화면에 표시되는 시각적 위치 설정
        # 불안정한 입자가 표시되지 않도록 하는 데 사용됨
        (self.visual_x_pos, self.visual_y_pos) = (self.x_pos, self.y_pos)

        # force 초기화
        (self.x_force, self.y_force) = (0.0, -G)

        # 오일러 적분(dt = 1)을 사용하여 속도 정의
        (self.x_vel, self.y_vel) = (
            self.x_pos - self.previous_x_pos,
            self.y_pos - self.previous_y_pos,
        )

        # 속도 계산
        velocity = sqrt(self.x_vel**2 + self.y_vel**2)

        # 속도가 너무 높으면 감소시킴
        if velocity > MAX_VEL:
            self.x_vel *= VEL_DAMP
            self.y_vel *= VEL_DAMP

        # 벽 제약 조건, 입자가 범위를 벗어나면 복원력을 생성하여 되돌림
        if self.x_pos < -SIM_W:
            self.x_force -= 0.3*(self.x_pos - -SIM_W) * WALL_DAMP
            self.visual_x_pos = -SIM_W

        # 댐 제약 조건, 댐이 dam에서 SIM_W로 이동하는 경우 적용, 사용 X
        if dam is True and self.x_pos > DAM:
            self.x_force -= (self.x_pos - DAM) * WALL_DAMP

        # 오른쪽 벽 제약 조건
        if self.x_pos > SIM_W:
            self.x_force -= 0.3*(self.x_pos - SIM_W) * WALL_DAMP
            self.visual_x_pos = SIM_W

        # 바닥 제약 조건
        if self.y_pos < BOTTOM:
            # 입자가 너무 낮아지지 않도록 BOTTOM 대신 SIM_W를 사용
            self.y_force -= 0.3*(self.y_pos - SIM_W) * WALL_DAMP
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