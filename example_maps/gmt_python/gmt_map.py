#!/usr/bin/env python3


import pygmt


fig = pygmt.Figure()
fig.coast(region=[190,206,68,73], shorelines=True, land="lightgreen", water="lightblue",projection="L200/65/60/77/12c",frame=["afg", "+tChukchi Moorings"])
fig.plot(x=-163.125, y=70.838, style="c0.3c", color="red", pen="black")
fig.plot(x=-164.223, y=71.231, style="c0.3c", color="red", pen="black")
fig.plot(x=-166.070, y=71.828, style="c0.3c", color="red", pen="black")
fig.plot(x=-160.514, y=71.038, style="c0.3c", color="red", pen="black")
fig.plot(x=-158.011, y=71.203, style="c0.3c", color="red", pen="black")
fig.text(text=["C1","C2","C3","C4","C5"], x=[-163.125, -164.223, -166.070, -160.514, -158.011],
         y=[71.038, 71.431, 72.028, 71.238, 71.403])
fig.savefig("chukchi_mooring_sites.png")
