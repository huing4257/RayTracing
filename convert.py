from PIL import Image

im = Image.open("../news.ppm")
im.save("../news.jpg")
