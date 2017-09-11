import os

json_text_template = """
{
     'empty_sphere_$f':
     {  environment:  { width:$f+1, height: 2*$f+1 }, // {pixmap: "empty_small.png"},
        plot_options: { plot_scale: 160/width,
                        line_thickness: 0.2,
                        vertex_size: 0.1,
                        save_image_interval: -1,
                        do_animation: true,
                        pause_for_visualization: false },
        print_options: { path_output_file: "outfiles/sphere_paths.data" },
        coord_scale: { chart_limits: "pixel_center", // "pixel_center" or "pixel_boundary"
                       x_min: 0.0, x_max: pi,        // x=latitude, y=longitude
                       y_min: 0.0, y_max: 2.0*pi },
        start:     [pi/2, pi],
        goal:      [pi/4, pi/2], // [5*pi/18, 
        cost:  { metric: [ "1.0",  "0.0",
                           "0.0",  "sin(x)*sin(x)" ], 
                 fun_to_pixmap: false },
     } ,
}
"""

algorithms = ["AStar", "ThetaStar", "SStar"]

for alg in algorithms:
    cmd = "cd ..; make nonuniform2d_PathPlanning DOSL_ALGORITHM="+alg
    print(cmd)
    os.system(cmd)
    for f in range(8, 81, 8):
        # create JSON file
        json_text = json_text_template.replace ("$f", str(f))
        json_file = open ("../files-expt/tmp_sphere.json", "w")
        json_file.write (json_text)
        json_file.close ()
        cmd = "../bin/nonuniform2d_PathPlanning ../files-expt/tmp_sphere.json empty_sphere_"+str(f)
        print(cmd)
        os.system(cmd)
        
