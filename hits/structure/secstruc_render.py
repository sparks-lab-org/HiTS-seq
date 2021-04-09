#https://gist.github.com/JoaoRodrigues/f9906b343d3acb38e39f2b982b02ecb0

import numpy as np
import matplotlib.patches as mpl_patches

def render_ss(ss_blocks, ax):
    # General settings
    helix_as_wave = False
    helix_as_cylinder = True

    base_color = 'grey'
    #base_color = 'firebrick'

    fc_sheet = base_color
    fc_helix = base_color
    fc_coil = base_color
    ec = 'none'

    width = 1.0 # General width controller
    helix_arc_width = 4.0
    coil_thickness_factor = 1/12.
    edge_thickness = 1.0
    sheet_thickness_factor = 2/3.
    turn_thickness_factor = 1/2.

    # Draw artists
    #width = 1.0
    for blk_idx, ss_blk in enumerate(ss_blocks):
        ss_type, start, last = ss_blk

        if helix_as_cylinder and ss_type == 'H':
            # Draw rectangle capped by two elipses
            # Elipse width = 1 residue
            # Origin is *center* of elipse

            # Order of drawing matters for overlap

            # First elipse
            elength = width / 2  # horizontal diameter
            height = width - 0.001  # vertical axis (slight offset bc of edge)
            origin = (start + elength/2, height/2)

            e = mpl_patches.Ellipse(origin, elength, height,
                                    linewidth=edge_thickness,
                                    edgecolor='black', facecolor=fc_helix)
            ax.add_patch(e)

            # Rectangle(s)
            length = last - start + 1 - elength  # deduct l of the ellipses
            height = width  # rectangle width is fraction of global
            origin = (start + elength/2, 0)  # origin is lower left: make it v-cntr

            e = mpl_patches.Rectangle(origin, length, height,
                                      edgecolor='none', facecolor=fc_helix)
            ax.add_patch(e)

            # Second elipse
            height = width - 0.001  # vertical axis
            origin = (last + 1 - elength/2, height/2)
            e = mpl_patches.Ellipse(origin, elength, height,
                                    linewidth=edge_thickness,
                                    edgecolor='black', facecolor=fc_helix)

            ax.add_patch(e)

        elif helix_as_wave and ss_type == 'H':
            # Draw as consecutive elyptical arcs
            height = width
            length = 0.5
            st_theta, en_theta = 0, 180
            for t_start in np.arange(start, last + 1, length):  # turns
                origin = (t_start + 0.25, height/2)
                e = mpl_patches.Arc(origin, length, height,
                                    linewidth=helix_arc_width,
                                    # Add a bit to each angle to avoid sharp cuts
                                    # that show as white lines in plot
                                    theta1=st_theta - 1, theta2=en_theta + 1,
                                    edgecolor=fc_helix)

                st_theta += 180
                en_theta += 180
                ax.add_patch(e)

        elif ss_type == 'E':
            # Draw arrow
            length = last - start + 1
            tail_height = width * sheet_thickness_factor
            head_height = width
            #print length, width/2.0, start, tail_height, width

            e = mpl_patches.FancyArrow(start, width/2.0,  # x, y of tail
                                       length, 0,  # dx, dy=0 -> flat arrow
                                       length_includes_head=True,
                                       head_length=length/4.0,
                                       head_width=head_height - 0.001,
                                       width=tail_height,
                                       facecolor=fc_sheet,
                                       edgecolor=ec,
                                       linewidth=edge_thickness)
            ax.add_patch(e)

        elif ss_type == 'T':
            # Draw turn as thin arc
            height = width
            length = last - start + 1
            st_theta, en_theta = 0, 180
            origin = (start + length / 2, height/2)
            e = mpl_patches.Arc(origin, length, height,
                                linewidth=helix_arc_width,
                                theta1=st_theta, theta2=en_theta,
                                edgecolor=fc_helix)

            ax.add_patch(e)

        else:  # draw line (thin Rectangle)
            length = last - start + 1
            height = width * coil_thickness_factor

            # Offset ends to fake continuity
            prev_blk_type = ss_blocks[blk_idx - 1][0]
            if prev_blk_type in ('H', 'T'):
                # Rougly the same size as the linewidth
                start -= 4/72
                length += 4/72
            elif prev_blk_type == 'E':
                # Go wild
                start -= 0.5
                length += 0.5

            if (blk_idx + 1) < len(ss_blocks):
                next_blk_type = ss_blocks[blk_idx + 1][0]
                if next_blk_type in ('H', 'T'):
                    length += 4/72
                elif next_blk_type == 'E':
                    length += 0.5

            origin = (start, width/2 - height/2)  # vertical center

            e = mpl_patches.Rectangle(origin, length, height,
                                      linewidth=edge_thickness,
                                      edgecolor=ec, facecolor=fc_coil)
            ax.add_patch(e)

        ax.set_yticklabels([])
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)
