import requests

from ipamd.public.utils.output import error
from ipamd.public.models.sequence import Sequence
from ipamd.public.models.data import Scalar
def func(seq: Sequence, **kwargs):
    """调用MolPhase的processSequence.php端点进行预测"""
    seq_string = seq.sequence
    try:
        response = requests.post(
            "https://molphase.sbs.ntu.edu.sg/processSequence.php",
            data={'sequence': seq_string},
            timeout=30
        )
        response.raise_for_status()  # 如果状态码不是200，抛出异常
        result = response.json()
        score = float(result['molphase_score'])
        return Scalar("MolPhase", score, "")
    except requests.exceptions.RequestException as e:
        error(f"请求出错: {e}")
        return None
