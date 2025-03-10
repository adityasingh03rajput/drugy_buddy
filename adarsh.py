import pygame
import socket
import pickle
import threading

# Game Setup
WIDTH, HEIGHT = 800, 600
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
FONT_SIZE = 24
SERVER_IP = '192.168.50.185'  # Change this to your IPv4 Address
SERVER_PORT = 5555
MAX_PLAYERS = 2

pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Soccer Shooter Game")
clock = pygame.time.Clock()
font = pygame.font.Font(None, FONT_SIZE)

# Game State
players = [{"x": 100, "y": 300, "color": (0, 0, 255)}, {"x": 700, "y": 300, "color": (255, 0, 0)}]
ball = {"x": 400, "y": 300, "velocity": [0, 0]}
scores = [0, 0]
game_state_lock = threading.Lock()

def draw_text(text, x, y, color=WHITE):
    """Draws text on the screen."""
    text_surface = font.render(text, True, color)
    screen.blit(text_surface, (x, y))

def run_server():
    """Handles the server logic."""
    server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server.bind((SERVER_IP, SERVER_PORT))
    server.listen(MAX_PLAYERS)
    print(f"Server started on {SERVER_IP}:{SERVER_PORT}, waiting for players...")

    connections = []
    
    def handle_client(conn, player_id):
        global ball, scores
        conn.send(pickle.dumps(player_id))  # Send player ID

        while True:
            try:
                data = conn.recv(1024)
                if not data:
                    break
                player_input = pickle.loads(data)

                with game_state_lock:
                    players[player_id]["x"] += player_input["dx"]
                    players[player_id]["y"] += player_input["dy"]

                    # Ball physics
                    ball["x"] += ball["velocity"][0]
                    ball["y"] += ball["velocity"][1]
                    ball["velocity"][0] *= 0.97
                    ball["velocity"][1] *= 0.97

                    # Goal detection
                    if ball["x"] < 0:
                        scores[1] += 1
                        ball = {"x": 400, "y": 300, "velocity": [0, 0]}
                    elif ball["x"] > 800:
                        scores[0] += 1
                        ball = {"x": 400, "y": 300, "velocity": [0, 0]}

                    game_state = {"players": players, "ball": ball, "scores": scores}
                    conn.send(pickle.dumps(game_state))

            except Exception as e:
                print(f"Player {player_id} disconnected: {e}")
                break

        conn.close()

    for i in range(MAX_PLAYERS):
        conn, addr = server.accept()
        print(f"Player {i+1} connected from {addr}")
        connections.append(conn)
        threading.Thread(target=handle_client, args=(conn, i)).start()

def run_client():
    """Handles the client logic."""
    conn = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    conn.connect((SERVER_IP, SERVER_PORT))
    player_id = pickle.loads(conn.recv(1024))
    print(f"Connected as Player {player_id + 1}")

    running = True
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        keys = pygame.key.get_pressed()
        movement = {"dx": 0, "dy": 0}

        if keys[pygame.K_LEFT]: movement["dx"] = -5
        if keys[pygame.K_RIGHT]: movement["dx"] = 5
        if keys[pygame.K_UP]: movement["dy"] = -5
        if keys[pygame.K_DOWN]: movement["dy"] = 5

        conn.send(pickle.dumps(movement))
        game_state = pickle.loads(conn.recv(1024))

        screen.fill(BLACK)
        for i, player in enumerate(game_state["players"]):
            pygame.draw.rect(screen, player["color"], (player["x"], player["y"], 40, 40))
        pygame.draw.circle(screen, WHITE, (game_state["ball"]["x"], game_state["ball"]["y"]), 15)
        draw_text(f"Score: {game_state['scores'][0]} - {game_state['scores'][1]}", 350, 20)

        pygame.display.flip()
        clock.tick(30)

if __name__ == "__main__":
    choice = input("Enter 's' to start server, 'c' to start client: ").strip().lower()
    if choice == 's':
        run_server()
    elif choice == 'c':
        run_client()
