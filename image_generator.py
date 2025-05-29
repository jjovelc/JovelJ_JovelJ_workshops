import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import random

def generate_landscape(width=1000, height=800, seed=None):
    """Generate a landscape illustration with mountains, trees, and sky"""
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(width/100, height/100), dpi=100)
    
    # Set background color (sky)
    time_of_day = random.choice(['day', 'sunset', 'night'])
    
    if time_of_day == 'day':
        sky_color = (0.53, 0.81, 0.98)  # Light blue
        sun_color = (1.0, 0.95, 0.1)    # Yellow
        mountain_colors = [(0.5, 0.55, 0.5), (0.4, 0.45, 0.4), (0.3, 0.35, 0.3)]  # Green-gray mountains
        tree_colors = [(0.0, 0.5, 0.0), (0.0, 0.4, 0.0), (0.0, 0.6, 0.0)]  # Green trees
        ground_color = (0.76, 0.69, 0.5)  # Tan ground
    elif time_of_day == 'sunset':
        sky_color = (0.98, 0.65, 0.3)  # Orange-red
        sun_color = (0.98, 0.4, 0.0)   # Deep orange
        mountain_colors = [(0.4, 0.3, 0.35), (0.3, 0.25, 0.3), (0.25, 0.2, 0.25)]  # Purplish mountains
        tree_colors = [(0.0, 0.3, 0.0), (0.0, 0.25, 0.0), (0.0, 0.4, 0.0)]  # Darker green trees
        ground_color = (0.65, 0.5, 0.4)  # Darker tan ground
    else:  # night
        sky_color = (0.05, 0.05, 0.2)  # Dark blue
        sun_color = (0.9, 0.9, 0.9)    # White/gray (moon)
        mountain_colors = [(0.2, 0.2, 0.25), (0.15, 0.15, 0.2), (0.1, 0.1, 0.15)]  # Dark blue-gray mountains
        tree_colors = [(0.0, 0.15, 0.0), (0.0, 0.1, 0.0), (0.0, 0.2, 0.0)]  # Very dark green trees
        ground_color = (0.2, 0.2, 0.3)  # Dark ground
    
    # Fill background
    ax.add_patch(patches.Rectangle((0, 0), width, height, facecolor=sky_color))
    
    # Add sun or moon
    sun_radius = random.uniform(40, 60)
    sun_x = random.uniform(sun_radius*1.5, width - sun_radius*1.5)
    sun_y = random.uniform(height * 0.6, height - sun_radius*1.5)
    sun = plt.Circle((sun_x, sun_y), sun_radius, color=sun_color)
    ax.add_patch(sun)
    
    # Add stars if night time
    if time_of_day == 'night':
        for _ in range(100):
            star_x = random.uniform(0, width)
            star_y = random.uniform(0.4 * height, height)
            star_size = random.uniform(1, 3)
            star_brightness = random.uniform(0.7, 1.0)
            star_color = (star_brightness, star_brightness, star_brightness)
            star = plt.Circle((star_x, star_y), star_size, color=star_color)
            ax.add_patch(star)
    
    # Add clouds
    if time_of_day != 'night':
        num_clouds = random.randint(3, 7)
        cloud_color = (0.95, 0.95, 0.95) if time_of_day == 'day' else (0.9, 0.7, 0.6)
        
        for _ in range(num_clouds):
            cloud_x = random.uniform(0, width)
            cloud_y = random.uniform(height * 0.6, height * 0.9)
            
            # Create a cloud with multiple overlapping circles
            num_puffs = random.randint(3, 6)
            base_radius = random.uniform(20, 40)
            
            for i in range(num_puffs):
                offset_x = random.uniform(-base_radius, base_radius)
                offset_y = random.uniform(-base_radius/2, base_radius/2)
                puff_radius = random.uniform(base_radius * 0.8, base_radius * 1.2)
                
                cloud_puff = plt.Circle((cloud_x + offset_x, cloud_y + offset_y), puff_radius, color=cloud_color)
                ax.add_patch(cloud_puff)
    
    # Generate mountains (3 layers for depth)
    mountain_heights = [0.5, 0.4, 0.3]  # Percentage of canvas height
    
    for i, height_factor in enumerate(mountain_heights):
        # Create a jagged mountain silhouette
        mountain_top = height * height_factor
        num_points = random.randint(15, 25)
        
        # Generate x coordinates
        x_coords = np.linspace(0, width, num_points)
        
        # Generate y coordinates (mountain peaks and valleys)
        y_coords = []
        for j in range(num_points):
            if j == 0 or j == num_points - 1:
                # Start and end at the bottom
                y_coords.append(0)
            else:
                # Random height for peaks and valleys
                variance = random.uniform(0.5, 1.0) if j % 2 == 0 else random.uniform(0.7, 0.9)
                y_coords.append(mountain_top * variance)
        
        # Create the mountain polygon
        mountain_points = np.column_stack([x_coords, y_coords])
        
        # Add bottom corners to close the polygon
        mountain_points = np.vstack([mountain_points, [width, 0], [0, 0]])
        
        # Create and add the mountain patch
        mountain = patches.Polygon(mountain_points, closed=True, facecolor=mountain_colors[i])
        ax.add_patch(mountain)
    
    # Add a ground/field
    ground_height = height * 0.2  # 20% of the canvas height
    ax.add_patch(patches.Rectangle((0, 0), width, ground_height, facecolor=ground_color))
    
    # Add trees
    num_trees = random.randint(15, 30)
    
    for _ in range(num_trees):
        # Determine tree position
        tree_x = random.uniform(0, width)
        tree_y = random.uniform(0, ground_height * 0.9)
        
        # Determine tree size
        tree_height = random.uniform(50, 120)
        trunk_width = tree_height * 0.1
        
        # Draw trunk
        trunk_color = (0.55, 0.27, 0.07)  # Brown
        trunk = patches.Rectangle(
            (tree_x - trunk_width/2, tree_y),
            trunk_width, tree_height * 0.4,
            facecolor=trunk_color
        )
        ax.add_patch(trunk)
        
        # Determine foliage type (pine or deciduous)
        tree_type = random.choice(['pine', 'deciduous'])
        
        if tree_type == 'pine':
            # Pine tree (triangular foliage)
            foliage_width = tree_height * 0.6
            foliage_color = random.choice(tree_colors)
            
            # Create 2-3 layers of triangles
            num_layers = random.randint(2, 3)
            layer_height = tree_height * 0.6 / num_layers
            
            for i in range(num_layers):
                layer_width = foliage_width * (1 - i * 0.2)
                layer_bottom = tree_y + tree_height * 0.3 + i * layer_height
                
                triangle_points = [
                    (tree_x, layer_bottom + layer_height),  # top
                    (tree_x - layer_width/2, layer_bottom),  # bottom left
                    (tree_x + layer_width/2, layer_bottom)   # bottom right
                ]
                
                foliage = patches.Polygon(triangle_points, closed=True, facecolor=foliage_color)
                ax.add_patch(foliage)
                
        else:
            # Deciduous tree (circular/cloud-like foliage)
            foliage_radius = tree_height * 0.3
            foliage_color = random.choice(tree_colors)
            foliage_center_y = tree_y + tree_height * 0.4 + foliage_radius * 0.8
            
            # Create main foliage circle
            main_foliage = plt.Circle((tree_x, foliage_center_y), foliage_radius, color=foliage_color)
            ax.add_patch(main_foliage)
            
            # Add some smaller circles for texture
            num_puffs = random.randint(3, 5)
            for i in range(num_puffs):
                puff_radius = foliage_radius * random.uniform(0.4, 0.7)
                angle = random.uniform(0, 2 * np.pi)
                distance = foliage_radius * random.uniform(0.5, 0.8)
                
                puff_x = tree_x + np.cos(angle) * distance
                puff_y = foliage_center_y + np.sin(angle) * distance
                
                puff = plt.Circle((puff_x, puff_y), puff_radius, color=foliage_color)
                ax.add_patch(puff)
    
    # Set axis limits and remove ticks
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.set_aspect('equal')
    ax.axis('off')
    
    return fig, ax

def save_landscape(filename="landscape.png", width=1000, height=800, seed=None):
    """Generate and save a landscape illustration"""
    fig, ax = generate_landscape(width, height, seed)
    plt.tight_layout(pad=0)
    plt.savefig(filename, bbox_inches='tight', pad_inches=0, dpi=100)
    plt.close(fig)
    print(f"Landscape saved as {filename}")

# Execute the script
if __name__ == "__main__":
    # You can specify a seed for reproducible results
    # or leave it as None for a random landscape each time
    random_seed = random.randint(1, 10000)
    print(f"Using seed: {random_seed}")
    save_landscape(filename="generated_landscape.png", seed=random_seed)
