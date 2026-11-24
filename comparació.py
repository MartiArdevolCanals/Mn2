import matplotlib.pyplot as plt
from unittest.mock import patch

with patch("matplotlib.pyplot.show"), patch("builtins.print"):
    from explícit import x
    from explícit import errors3
    from implicit import errors1 as eimp
    from crank_nicolson import errors2

plt.close("all")

plt.figure()
plt.plot(x,errors3, label="Error per Euler explícit $\gamma$=0.25")
plt.plot(x,eimp, label="Error per Euler implíct $\gamma$=0.5")
plt.plot(x,errors2, label="Error per Crank-Nicolson $\gamma$=0.5")
plt.xlabel("Posició (cm)")
plt.ylabel("Error absolut (°C)")
plt.title("Comparació dels errors dels mètodes numèrics")
plt.legend()
plt.show()