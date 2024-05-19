import uuid
import os

def generate_id_file(name, resolution, folder):
    while True:
        unique_id = str(uuid.uuid4().hex[:10])
        output_file_id = f"{name}_{unique_id}.{resolution}"
        output_file_id = folder + str(output_file_id)
        print('try')

        if not os.path.exists(output_file_id):
            return output_file_id