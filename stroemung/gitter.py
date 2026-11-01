#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Das versetzte Gitter (MAC) -- die topologisch richtige Aufstellung.

Im Netzwerkbild leben POTENTIALE auf Knoten und FLUESSE auf Zweigen. Der Druck ist ein
Potential, der Durchfluss ein Fluss. Also:

    p   in der Zellmitte          Knoten      (nx, ny)
    u   auf der senkrechten Flaeche   Zweig   (nx+1, ny)
    v   auf der waagerechten Flaeche  Zweig   (nx, ny+1)

Das ist nicht die Reparatur eines Fehlers, sondern die richtige Zuordnung. Legt man die
Geschwindigkeit auf denselben Knoten wie den Druck, wird der Druckgradient zentral gebildet
und ein Schachbrettfeld p_i = (-1)^i A ist unsichtbar:

    zentral:   (p_{i+1} - p_{i-1}) / 2h = 0
    versetzt:  (p_{i+1} - p_i) / h     = -2A/h     <- sofort bestraft

Harlow und Welch 1965.
"""
FASSUNG = "1.0"
import numpy as np


class Gitter:
    """Versetztes Gitter auf [0,Lx] x [0,Ly] mit nx x ny Zellen."""

    def __init__(self, Lx, Ly, nx, ny):
        self.Lx, self.Ly, self.nx, self.ny = float(Lx), float(Ly), int(nx), int(ny)
        self.dx = self.Lx / self.nx
        self.dy = self.Ly / self.ny
        # Zellmitten (Druckknoten)
        self.xp = (np.arange(nx) + 0.5) * self.dx
        self.yp = (np.arange(ny) + 0.5) * self.dy
        # Flaechenmitten (Geschwindigkeitszweige)
        self.xu = np.arange(nx + 1) * self.dx          # u sitzt bei x = i*dx
        self.yu = self.yp                              # ... und in Zellmitte in y
        self.xv = self.xp
        self.yv = np.arange(ny + 1) * self.dy          # v sitzt bei y = j*dy

    def leer_u(self):
        return np.zeros((self.nx + 1, self.ny))

    def leer_v(self):
        return np.zeros((self.nx, self.ny + 1))

    def leer_p(self):
        return np.zeros((self.nx, self.ny))

    def divergenz(self, u, v):
        """Knotenregel je Zelle: Summe der Zweigstroeme durch die vier Flaechen.

        Geteilt durch das Zellvolumen, also die Divergenz. Null heisst: Kirchhoff erfuellt.
        """
        return ((u[1:, :] - u[:-1, :]) / self.dx
                + (v[:, 1:] - v[:, :-1]) / self.dy)

    def __repr__(self):
        return (f"Gitter {self.nx}x{self.ny} auf {self.Lx}x{self.Ly}, "
                f"dx={self.dx:.5g}, dy={self.dy:.5g}")


if __name__ == "__main__":
    g = Gitter(1.0, 0.5, 8, 4)
    print(g)
    print(f"  p {g.leer_p().shape}   u {g.leer_u().shape}   v {g.leer_v().shape}")
    # Probe: ein divergenzfreies Feld muss Divergenz null geben.
    # u = x, v = -y ist divergenzfrei (du/dx + dv/dy = 1 - 1 = 0)
    X, Y = np.meshgrid(g.xu, g.yu, indexing="ij")
    u = X.copy()
    X2, Y2 = np.meshgrid(g.xv, g.yv, indexing="ij")
    v = -Y2.copy()
    d = g.divergenz(u, v)
    print(f"  Probe u=x, v=-y (divergenzfrei): groesste Divergenz {np.abs(d).max():.3e}")
    print(f"  [{'RICHTIG' if np.abs(d).max() < 1e-13 else 'FALSCH'}] Knotenregel stimmt")
