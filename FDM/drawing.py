import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
from tkinter import *

class GraphDrawer:
  def __init__(self):
    self.fig, self.ax = plt.subplots()
    self.points = {'x': [], 'y': []}
    self.is_drawing = False  # Flag to track drawing state

    self.fig.canvas.mpl_connect('button_press_event', self.on_click)
    self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
    self.fig.canvas.mpl_connect('button_release_event', self.on_release)

  def on_click(self, event):
    if event.button == MouseButton.LEFT:
      self.is_drawing = True  # Start drawing on left button press
      self.points['x'].append(event.xdata)
      self.points['y'].append(event.ydata)
      self.draw()

  def on_motion(self, event):
    if self.is_drawing and event.button == MouseButton.LEFT:
      # Update last point while left button is pressed and dragged
      self.points['x'][-1] = event.xdata
      self.points['y'][-1] = event.ydata
      self.draw()

  def on_release(self, event):
    if event.button == MouseButton.LEFT:
      self.is_drawing = False  # Stop drawing on left button release
      self.draw()

  def draw(self):
    self.ax.clear()
    self.ax.plot(self.points['x'], self.points['y'], marker='o')
    self.ax.set_title('Draw a Graph')
    self.ax.set_xlabel('X')
    self.ax.set_ylabel('Y')
    self.fig.canvas.draw()

  def show(self):
    plt.show()

def main():
  root = Tk()
  root.withdraw()  # Hide the default tkinter window

  graph_drawer = GraphDrawer()
  graph_drawer.show()

if __name__ == "__main__":
  main()
