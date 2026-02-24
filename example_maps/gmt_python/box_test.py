#!/usr/bin/env python3


import pygmt


fig = pygmt.Figure()
fig.coast(region=[-169.5,-168.5,65.5,66], shorelines=True, land="lightgreen", water="lightblue",projection="L200/65/60/77/12c",frame=["afg", "+tChukchi Moorings"])
fig.plot(x=-169, y=65.75, style="c0.3c", color="red", pen="black")
fig.plot(x=[-169.1,168.9], y=[65.65, 65.85], style="c0.3c", color="red", pen="black")
fig.savefig("box_test.png")
