import pygame
import socket
import threading
import pickle
import random

# Constants
WIDTH, HEIGHT = 800, 600
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
FONT_SIZE = 24
SERVER_IP = '0.0.0.0'  # Host listens on all interfaces
SERVER_PORT = 5555

# Initialize Pygame
pygame.init()
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Soccer Shooter Game")
clock = pygame.time.Clock()
font = pygame.font.Font(None, FONT_SIZE)

# Game state
players = []
ball = {"x": WIDTH // 2, "y": HEIGHT // 2, "velocity": [0, 0]}
scores = [0, 0]
game_state_lock = threading.Lock()


def generate_room_code():
    """Generate a random 4-digit room code."""
    return str(random.randint(1000, 9999))


def start_server(room_code):
    """Start the server and wait for clients to join."""
    server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server.bind((SERVER_IP, SERVER_PORT))
    server.listen(2)
    print(f"Room code: {room_code}. Waiting for players to join...")

    conn1, addr1 = server.accept()
    print(f"Player 1 connected from {addr1}")
    conn2, addr2 = server.accept()
    print(f"Player 2 connected from {addr2}")
    return conn1, conn2


def handle_client(conn, player_id):
    """Handle communication with a client."""
    while True:
        try:
            # Receive player input from the client
            data = conn.recv(1024)
            if not data:
                break
            player_input = pickle.loads(data)

            # Update game state based on player input
            with game_state_lock:
                players[player_id]["x"] += player_input["dx"]
                players[player_id]["y"] += player_input["dy"]
                if player_input["shoot"]:
                    # Handle shooting logic (e.g., create bullets)
                    pass

                # Update ball physics
                ball["x"] += ball["velocity"][0]
                ball["y"] += ball["velocity"][1]
                ball["velocity"][0] *= 0.97
                ball["velocity"][1] *= 0.97

                # Check for goals
                if ball["x"] < 0:
                    scores[1] += 1
                    ball["x"], ball["y"] = WIDTH // 2, HEIGHT // 2
                    ball["velocity"] = [0, 0]
                elif ball["x"] > WIDTH:
                    scores[0] += 1
                    ball["x"], ball["y"] = WIDTH // 2, HEIGHT // 2
                    ball["velocity"] = [0, 0]

                # Send updated game state back to the client
                game_state = {
                    "players": players,
                    "ball": ball,
                    "scores": scores
                }
                conn.send(pickle.dumps(game_state))
        except Exception as e:
            print(f"Error: {e}")
            break
    conn.close()


def draw_text(text, x, y, color=WHITE):
    """Draw text on the screen."""
    text_surface = font.render(text, True, color)
    screen.blit(text_surface, (x, y))


def main():
    # Game setup
    role = None
    room_code = None
    conn = None

    # Main menu
    while role is None:
        screen.fill(BLACK)
        draw_text("1. Host Party", WIDTH // 2 - 50, HEIGHT // 2 - 50)
        draw_text("2. Join Party", WIDTH // 2 - 50, HEIGHT // 2)
        pygame.display.flip()

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                return
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_1:
                    role = "host"
                elif event.key == pygame.K_2:
                    role = "join"

    # Host Party
    if role == "host":
        room_code = generate_room_code()
        conn1, conn2 = start_server(room_code)
        players.append({"x": 100, "y": HEIGHT // 2, "color": (0, 0, 255), "health": 100})
        players.append({"x": WIDTH - 100, "y": HEIGHT // 2, "color": (255, 0, 0), "health": 100})

        # Start threads to handle clients
        thread1 = threading.Thread(target=handle_client, args=(conn1, 0))
        thread2 = threading.Thread(target=handle_client, args=(conn2, 1))
        thread1.start()
        thread2.start()
        thread1.join()
        thread2.join()

    # Join Party
    elif role == "join":
        screen.fill(BLACK)
        draw_text("Enter Room Code:", WIDTH // 2 - 100, HEIGHT // 2 - 50)
        input_text = ""
        while True:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    pygame.quit()
                    return
                if event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_RETURN:
                        room_code = input_text
                        break
                    elif event.key == pygame.K_BACKSPACE:
                        input_text = input_text[:-1]
                    else:
                        input_text += event.unicode

            screen.fill(BLACK)
            draw_text("Enter Room Code:", WIDTH // 2 - 100, HEIGHT // 2 - 50)
            draw_text(input_text, WIDTH // 2 - 50, HEIGHT // 2)
            pygame.display.flip()

            if room_code:
                break

        # Connect to server
        conn = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        conn.connect((SERVER_IP, SERVER_PORT))
        player_id = int(conn.recv(1024).decode())  # Receive player ID from the server
        print(f"You are Player {player_id + 1}")

        # Game loop
        running = True
        while running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False

            # Send player input to the server
            keys = pygame.key.get_pressed()
            dx, dy = 0, 0
            shoot = False
            if keys[pygame.K_LEFT]:
                dx -= 3
            if keys[pygame.K_RIGHT]:
                dx += 3
            if keys[pygame.K_UP]:
                dy -= 3
            if keys[pygame.K_DOWN]:
                dy += 3
            if keys[pygame.K_SPACE]:
                shoot = True

            player_input = {"player_id": player_id, "dx": dx, "dy": dy, "shoot": shoot}
            conn.send(pickle.dumps(player_input))

            # Receive updated game state from the server
            data = conn.recv(4096)
            game_state = pickle.loads(data)

            # Render the game state
            screen.fill(BLACK)
            for player in game_state["players"]:
                pygame.draw.circle(screen, player["color"], (player["x"], player["y"]), 20)
            pygame.draw.circle(screen, (255, 255, 0), (int(game_state["ball"]["x"]), int(game_state["ball"]["y"])), 15)
            draw_text(f"Score: {game_state['scores'][0]} - {game_state['scores'][1]}", WIDTH // 2 - 50, 50)
            pygame.display.flip()
            clock.tick(60)


# Ensure the script runs properly
if __name__ == "__main__":
    main()
