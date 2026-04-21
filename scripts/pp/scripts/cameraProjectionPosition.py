import numpy as np

# Camera parameters
camera_pos = np.array([-31.5952148608975, 24.9751367723974, 129.4358263432127])
focal_point = np.array([235.96773270122878, -15.453363381193055, -183.05801613248408])
view_up = np.array([0.03300926844990681, 0.9944008217976584, -0.10038622318099039])

# Additional camera parameters
view_angle = 30  # Field of view in degrees
eye_angle = 2    # Stereo eye angle (not typically used for orientation axes)
parallel_scale = 1.0  # For parallel projection (if applicable)

# Normalize view_up vector
view_up = view_up / np.linalg.norm(view_up)

# Compute camera coordinate system
view_dir = focal_point - camera_pos
distance = np.linalg.norm(view_dir)
view_dir = view_dir / distance  # normalize

# Right vector (camera's x-axis in world space)
right = np.cross(view_dir, view_up)
right = right / np.linalg.norm(right)

# Up vector (camera's y-axis in world space, recomputed for orthogonality)
up = np.cross(right, view_dir)
up = up / np.linalg.norm(up)

print("Camera basis vectors:")
print(f"View direction: {view_dir}")
print(f"Right: {right}")
print(f"Up: {up}")
print(f"Distance to focal point: {distance:.2f}")

# 3D world axes
x_axis = np.array([1, 0, 0])
y_axis = np.array([0, 1, 0])
z_axis = np.array([0, 0, 1])

# Project 3D vector onto 2D screen plane
def project_to_2d(vector_3d, scale=1.0):
    """
    Project a 3D world vector onto the 2D camera screen.
    Returns (horizontal, vertical) coordinates in screen space.
    """
    # For orientation axes, we typically use orthographic projection
    # of the direction vectors themselves (not positioned in 3D space)
    x_2d = np.dot(vector_3d, right) * scale
    y_2d = np.dot(vector_3d, up) * scale
    return x_2d, y_2d

# Alternative: Perspective projection (if needed)
def project_perspective(point_3d, scale=1.0):
    """
    Full perspective projection of a 3D point.
    """
    # Transform to camera coordinates
    rel_pos = point_3d - camera_pos
    
    # Project onto camera plane
    cam_x = np.dot(rel_pos, right)
    cam_y = np.dot(rel_pos, up)
    cam_z = np.dot(rel_pos, view_dir)  # depth
    
    # Perspective division (if cam_z > 0)
    if cam_z > 0.001:
        # Field of view scaling factor
        fov_rad = np.deg2rad(view_angle)
        aspect = 1.0  # assume square viewport
        
        # Perspective projection
        screen_x = (cam_x / cam_z) * scale
        screen_y = (cam_y / cam_z) * scale
        return screen_x, screen_y
    else:
        return 0, 0

# For orientation axes: project unit vectors with orthographic projection
# (This is standard for orientation widgets)
scale = 1.0  # Adjust for arrow length in TikZ

print("\n=== Orthographic Projection (Standard for Orientation Axes) ===")
x_proj = project_to_2d(x_axis, scale)
y_proj = project_to_2d(y_axis, scale)
z_proj = project_to_2d(z_axis, scale)

print(f"X-axis projection: ({x_proj[0]:+.3f}, {x_proj[1]:+.3f})")
print(f"Y-axis projection: ({y_proj[0]:+.3f}, {y_proj[1]:+.3f})")
print(f"Z-axis projection: ({z_proj[0]:+.3f}, {z_proj[1]:+.3f})")

print("\n=== TikZ Code ===")
print("\\begin{scope}[xshift=-4.95cm, yshift=-1.2cm]")
print(f"    % X-axis: red")
print(f"    \\draw[-{{Stealth[length=2mm]}}, red!80, line width=1.5pt]")
print(f"        (0,0) -- ({x_proj[0]:.3f}, {x_proj[1]:.3f}) node[right] {{$x$}};")
print(f"    ")
print(f"    % Y-axis: green")
print(f"    \\draw[-{{Stealth[length=2mm]}}, green!60!black, line width=1.5pt]")
print(f"        (0,0) -- ({y_proj[0]:.3f}, {y_proj[1]:.3f}) node[above left] {{$y$}};")
print(f"    ")
print(f"    % Z-axis: blue")
print(f"    \\draw[-{{Stealth[length=2mm]}}, blue!80, line width=1.5pt]")
print(f"        (0,0) -- ({z_proj[0]:.3f}, {z_proj[1]:.3f}) node[below right] {{$z$}};")
print("\\end{scope}")

# Optional: Show angles for verification
print("\n=== Verification ===")
angle_xy = np.rad2deg(np.arctan2(x_proj[1], x_proj[0]))
angle_yy = np.rad2deg(np.arctan2(y_proj[1], y_proj[0]))
angle_zy = np.rad2deg(np.arctan2(z_proj[1], z_proj[0]))
print(f"X-axis angle from horizontal: {angle_xy:.1f}°")
print(f"Y-axis angle from horizontal: {angle_yy:.1f}°")
print(f"Z-axis angle from horizontal: {angle_zy:.1f}°")
