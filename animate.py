import imageio
import os
from tqdm import tqdm
from typing import List, Union


def animate(input_images: Union[List[str], str], output: str, fps: int = 1, quiet: bool = False) -> str:
    """
    Creates an MP4 from still images.
    :param input_images: Either a list of paths to the images, or the path to the folder containing the images.  In the
    latter case, the images in the folder are sorted alphabetically.
    :param output: The location that the MP4 will be saved to.
    :param fps: The framerate of the MP4 in frames per second.  Default 1.
    :param quiet: Whether to suppress printing of progress.  Default False.
    :return: The path to the MP4.
    """

    if type(input_images) is str:
        input_images = sorted(os.listdir(input_images))

    writer = imageio.get_writer(output, fps=fps)

    if quiet:
        for im in input_images:
            writer.append_data(imageio.imread(im))
        writer.close()

    else:
        for im in tqdm(input_images, desc="Writing frames"):
            writer.append_data(imageio.imread(im))
        print("Closing writer...")
        writer.close()
        print("Finished animating.")

    return output
